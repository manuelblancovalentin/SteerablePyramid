function ThisReconstruction = SteerableReconstruction(ThisPyramid, Options, flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reconstructs the original image using all bands from all
% scales of the Steerable Pyramid in ThisPyramid. 
%
% Inputs:  ThisPyramid - Structure containing the Steerable Pyramid of the
%                        analyzed texture. This is the first output of the
%                        function 'SteerablePyramid'. In order to get it,
%                        use: 
%                           [ThisPyramid,TexFeats] = SteerablePyramid(ThisOriginalImage,Options)
%
%              Options - Structure with options for the code to work. These
%                        options are:
%
%                                   I) WindowSize -> Size of the Window
%                                                    that will be used for 
%                                                    extracting the
%                                                    textural features of
%                                                    the steerable pyramid.
%
%                                   II) NLevels -> Number of levels of the
%                                                 steerable pyramid. An
%                                                 integer number greater
%                                                 than 0 and smaller than 
%                                                 floor(log2(min(size(ThisOriginalImage))/WindowSize))
%
%                                   III) NBands -> Number of bands for each
%                                                  level of the Steerable
%                                                  Pyramid. Each band
%                                                  corresponds to a certain
%                                                  orientation. It must be 
%                                                  an integer value between
%                                                  1 and Inf.
%
%
%                 flag - Boolean value. If true, it displays the different
%                        processes progress in the display. If false, it
%                        remains silent during all the process.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Manuel Blanco Valentin
% Industrial Engineer & Project Analyst
% E-mail: mbvalentin@cbpf.br
%
% Centro Brasileiro de Pesquisas (CBPF) - CENPES - PETROBRAS
% Rio de Janeiro - Brazil - 2017
%
% LICENSE: This code is open-source, feel free to use it under your own
% responsability. Feel free to share this code, but please, do not delete
% these comments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin < 3, flag = true; end

if ~islogical(flag), flag = true; end;

res = 256;
twidth = 1;
X = (pi/(2*res)).*(-(res+1):1);
Y = cos(X).^2;
X = (2*twidth/pi).*X;

% We have to make a lil correction so the mask will have the values we
% want:
Y(1) = Y(2);
Y(end) = Y(end-1);

% In the original paper nBands --> 'K'
alphaK = 2^(Options.NBands-1)*factorial(Options.NBands-1)/sqrt(Options.NBands*factorial(2*(Options.NBands-1)));

% In the original code (Eero & Portilla), they call 'alfa' the angle:
% (theta - pi*k/(nBands-1)), in their paper (page 55 or 7, depending on the
% edition), which is very confusing actually, because they use alphaK as a 
% constant. Here theta is theta (in their equations) and the parameter that
% they call alfa is here called phi.
len = 1024;
theta = pi*(-(2*len+1):(len+1))./len;
% phi = mod(theta+pi,2*pi)-pi;

% In the original code they call Gk(theta) Ycosn. ALSO, THEY ADDED A FACTOR
% OF 2 FROM THEIR ORIGINAL EQUATION, IDK WHY.
%Gk = 2*alphaK*(cos(phi).^(Options.NBands-1)).*(abs(phi) < pi/2);

angleX = theta; 
angleY = alphaK.*cos(angleX).^(Options.NBands - 1);

ThisReconstruction = cellfun(@(x) zeros(size(x,1),size(x,2),1),ThisPyramid.LowPass,'un',0);

% add leftovers 
ThisReconstruction{end} = ThisPyramid.LowPass{end};

mayIprint(sprintf('\tReconstructing Steerable Pyramid:\n'),flag); 
for n=Options.NLevels:-1:1
    
    mayIprint(sprintf('\t\tLevel: %i\n',n),flag); tic;
    ThisSubPyr = ThisPyramid.LowPass{n};
    
    [M,N,~] = size(ThisSubPyr);

    %% 1) Basic Masks (lowpass, high pass...)
    
    centerj = floor(N/2)+1;
    centeri = floor(M/2)+1;

    % Meshgrid for the filter
    jj = ((1:N)-centerj)/(N/2);
    ii = ((1:M)-centeri)/(M/2);
    [xgrid,ygrid] = meshgrid(jj,ii);

    % extract log-distance and angle from the grid
    omega = atan2(ygrid,xgrid);
    rho = hypot(xgrid,ygrid);
    rho(centeri,centerj) = rho(centeri,centerj-1);
    rho = log2(rho);


    
    % 1.2) High pass mask (H0Mask) and high-pass component (H0). This component
    % will not be used in the building of the steerable pyramid.
    H0Mask_band = reshape(interp1(X-log2(2),sqrt(Y),rho(:), 'linear', 'extrap'),size(rho));

    Reconstruction_tmp = ThisReconstruction{n};
    for nb=1:Options.NBands
        mayIprint(sprintf('\t\t\tBand %i...\n',nb),flag);
        % angle mask
        anglemask = reshape(interp1(angleX(1)+pi*(nb-1)/Options.NBands + (angleX(2)-angleX(1))*(0:numel(angleY)-1),...
                         angleY,omega(:),'linear','extrap'),size(omega));
        Reconstruction_tmp = Reconstruction_tmp + ((1i)^(Options.NBands-1))*fftshift(fft2(real(ThisSubPyr(:,:,nb)))).*anglemask.*H0Mask_band;
    end

    % 1.1) Low pass mask (L0Mask) and Component (L0) we will use this low pass 
    % component L0 to build the steerable pyramid (check out the schema 
    % proposed by Portilla & Simoncelli in their article)
    L0Mask = reshape(interp1(X,sqrt(abs(1-Y)),rho(:), 'linear', 'extrap'),size(rho));
    Reconstruction_tmp = Reconstruction_tmp.*L0Mask;
    
    % 1.2) High pass mask (H0Mask) and high-pass component (H0). This component
    % will not be used in the building of the steerable pyramid.
    %H0Mask = reshape(interp1(X,sqrt(Y),rho(:), 'linear', 'extrap'),size(rho));
    %warning('this part is missing!')
    
    ThisReconstruction{n} = real(ifft2(ifftshift(Reconstruction_tmp))) + expand(ThisReconstruction{n+1},2)/4; % Current reconstruction plus previous reconstruction (notice that in the first iteration our previous reconstruction is actually the left-overs of the steering process, a.k.a. the LP left-overs).

    mayIprint(sprintf('\t\tProcess completed in: %f sec\n',toc),flag);
    
     % Drop X one octave
     %X = X + log2(2);
    
end


end

function mayIprint(string,flag)
if flag, fprintf('%s',string); end
end