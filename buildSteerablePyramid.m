% Extraction and Building of the Steerable Pyramid
function ThisPyramid = buildSteerablePyramid(ThisImage, Options, silent_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function build the Steerable pyramid with the options specified in
% the Options variable (Number of Levels/Scales, Number of Subbands, etc.)
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
%          silent_flag - Boolean value. If true, it displays the different
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


if nargin < 3, silent_flag = false; end
if ~islogical(silent_flag), silent_flag = false; end

[M,N,~] = size(ThisImage);


%% 1) Basic Masks (lowpass, high pass...)
res = 256;
twidth = 1;
X = (pi/(2*res)).*(-(res+1):1);
Y = cos(X).^2;
X = (2*twidth/pi).*X;

% We have to make a lil correction so the mask will have the values we
% want:
Y(1) = Y(2);
Y(end) = Y(end-1);

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


% 1.1) Low pass mask (L0Mask) and Component (L0) we will use this low pass 
% component L0 to build the steerable pyramid (check out the schema 
% proposed by Portilla & Simoncelli in their article)
L0Mask = reshape(interp1(X,sqrt(1-Y),rho(:), 'linear', 'extrap'),size(rho));
L0 = fftshift(fft2(ThisImage)).*L0Mask;

% 1.2) High pass mask (H0Mask) and high-pass component (H0). This component
% will not be used in the building of the steerable pyramid.
H0Mask = reshape(interp1(X,sqrt(Y),rho(:), 'linear', 'extrap'),size(rho));
H0 = fftshift(fft2(ThisImage)).*H0Mask;
H0 = ifft2(ifftshift(H0));


%% 2) Using the Low Pass component of the image (L0), let's create the 
% pyramid. The methodology should be, in theory, simple: 
%   a) Take the Low-pass component of the filtered image and use a series 
%      of oriented subbands (Bk) to filter the image. These components will
%      carry the information about orientation, gradient and direction of
%      the image textures.
%   b) Filter the initial low-pass component with another low-pass filter.
%   This new low-pass component will be the starting point of the next
%   step.

% In the original paper nBands --> 'K'
alphaK = 2^(Options.NBands-1)*factorial(Options.NBands-1)/sqrt(Options.NBands*factorial(2*(Options.NBands-1)));

% In the original code (Eero & Portilla), they call 'alfa' the angle:
% (theta - pi*k/(nBands-1)), in their paper (page 55 or 7, depending on the
% edition), which is very confusing actually, because they use alphaK as a 
% constant. Here theta is theta (in their equations) and the parameter that
% they call alfa is here called phi.
len = 1024;
theta = pi*(-(2*len+1):(len+1))./len;
phi = mod(theta+pi,2*pi)-pi;

% In the original code they call Gk(theta) Ycosn. ALSO, THEY ADDED A FACTOR
% OF 2 FROM THEIR ORIGINAL EQUATION, IDK WHY.
Gk = 2*alphaK*(cos(phi).^(Options.NBands-1)).*(abs(phi) < pi/2);

Lk = L0;

ThisPyramid = cell(1,Options.NLevels);
TheseLP = cell(1,Options.NLevels); % These are the LP marginal responses that we will use later when extracting some params

for n=1:Options.NLevels
    mayIprint(sprintf('\tBuilding Steerable Pyramid: Level %i\n',n),silent_flag); tic;
    
    % Downscale X by one octave
    X = X - log2(2);
    Hk = reshape(interp1(X(1) + (X(2)-X(1))*(0:numel(Y)-1),...
                         sqrt(Y),rho(:),'linear','extrap'),size(rho));
    
    subBands_tmp = [];
    for nb=1:Options.NBands
        mayIprint(sprintf('\t\tBand %i...\n',nb),silent_flag);
        % Creation of the subband filter (for a certain orientation)
        Bk = reshape(interp1(theta(1)+pi*(nb-1)/Options.NBands + (theta(2)-theta(1))*(0:numel(Gk)-1),...
                         Gk,omega(:),'linear','extrap'),size(omega));
        % Thus, the component Bk (B0, B1, B2...)
        Bk = ((-1i)^(Options.NBands-1)).*Lk.*Bk.*Hk;
        % Take inverse fourier to get the filtered image
        subBands_tmp = cat(3,subBands_tmp,ifft2(ifftshift(Bk)));
    end
    mayIprint(sprintf('\t\tComplete in (%f sec)\n',toc),silent_flag);
    
    ThisPyramid{n} = subBands_tmp;
    TheseLP{n} = ifft2(ifftshift(Lk));
    
    % Now let's downsample the image for the next step
    [M,N,~] = size(Lk);
    ist = ceil(M/2)-ceil(M/4)+1;
    ien = ceil(M/2)+ceil(M/4);
    jst = ceil(N/2)-ceil(N/4)+1;
    jen = ceil(N/2)+ceil(N/4);
   
    Lk = Lk(ist:ien,jst:jen,:);
    omega = omega(ist:ien,jst:jen,:);
    rho = rho(ist:ien,jst:jen,:);
    
    % create new low pass mask
    LkMask = reshape(interp1(X,abs(sqrt(1-sqrt(Y).^2)),rho(:), 'linear', 'extrap'),size(rho));
    Lk = Lk.*LkMask;
    
end

% Add NLevels + 1 filtering
TheseLP{end+1} = ifft2(ifftshift(Lk));
ThisPyramid{end+1} = real(TheseLP{end});

% Subtract mean from last low pass
ThisPyramid{end} = real(ThisPyramid{end}) - mean(real(ThisPyramid{end}(:)));

% Add the High-pass part of the response to the pyramid
ThisPyramid = struct('HighPass',{real(H0)},'LowPass',{ThisPyramid},'Marginals',{TheseLP});

% % We can display all subbands in a single image using the following
% % function
% dispPyramid(ThisImage,ThisPyramid.LowPass(1:end-1),ThisPyramid.Marginals(1:end-1))


end

function mayIprint(string,flag)
if flag, fprintf('%s',string); end
end