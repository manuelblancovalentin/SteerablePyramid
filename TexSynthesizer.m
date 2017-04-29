function [ThisSynthesized, ThisEvolution] = TexSynthesizer(ThisPyramid,TexFeats,Options,ThisMask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a synthetic image based on the Pyramid
% (ThisPyramid) extracted using the function 'SteerablePyramid', the
% textural features extracted from each level and subband of this Pyramid
% (TexFeats) and a series of Options (Size of the output image and Number
% of Iterations). It is also possible to specify a Mask (ThisMask), which
% is an optional argument, so that a certain part of the original image is
% always imposed to the synthesized image at each iteration. This last
% process produces an extrapolation of the original texture (kind of
% 'filling' the blank spaces specified with zeros at ThisMask).
%
% Inputs:  ThisPyramid - Structure containing the Steerable Pyramid of the
%                        analyzed texture. This is the first output of the
%                        function 'SteerablePyramid'. In order to get it,
%                        use: 
%                           [ThisPyramid,TexFeats] = SteerablePyramid(ThisOriginalImage,Options)
% 
%             TexFeats - Structure containing the Textural Feats obtained
%                        from all levels and scales of the Steerable 
%                        Pyramid. This is the second output of the function
%                        'Steerable Pyramid'. These features are essential
%                        to create the synthetic image, as we will impose
%                        them later to recreate the textural properties of
%                        the original image, and make the synthetic image
%                        match them. In order to get it, use the same line
%                        of code written in the previous paragraph (the
%                        same used for obtaining ThisPyramid).
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
%                                   IV) Size -> Size of the output
%                                               synthetic image specified
%                                               as a vector, e.g.: [100,200]
%
%
%                                   V) MaxIter -> Maximum Number of
%                                                 iterations for the
%                                                 synthesis process. Must
%                                                 be an integer between 1
%                                                 and Inf.
%
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



if isempty(ThisPyramid) || ~isstruct(ThisPyramid)
    error('Empty Pyramid.');
end

% check options
Options = checkInputOptions(Options);

% Variable definitions
WindowSize = ThisPyramid.Options.WindowSize;
NLevels = ThisPyramid.Options.NLevels;
NBands = ThisPyramid.Options.NBands;

% function definitions
dPhase = @(M) abs(M).*exp( (2i).*atan2(real(M),imag(M) ) );

%% 1) Create a gaussian white noise with the desired size. This noise will 
% be the seed of our Synthesized image. The initialized noise will have the
% same mean and variance of the original image (mu0, sigma0) 
mu0 = TexFeats.MarginalStatistics.Original(1);
sigma0 = TexFeats.MarginalStatistics.Original(2);
skew0 = TexFeats.MarginalStatistics.Original(3);
kurt0 = TexFeats.MarginalStatistics.Original(4);
min0 = TexFeats.MarginalStatistics.Original(5);
max0 = TexFeats.MarginalStatistics.Original(6);
% Mean of Magnitudes
magMeans0 = TexFeats.MarginalStatistics.MagnitudeMeans;
% Variance of high-frequency band
HPvar0 = TexFeats.MarginalStatistics.HP;
% Skewness and kurtosis of the lowest LP band
skewLP = TexFeats.MarginalStatistics.LP.Skewness;
kurtLP = TexFeats.MarginalStatistics.LP.Kurtosis;

% RawCoefficientCorrelation (in the original code called autoCorrReal),
% which we will impose to the synthesized image
rawCorr0 = TexFeats.RawCoefficientCorrelation;
% Auto Correlation of Magnitude Responses
autoCrossMag0 = TexFeats.CoefficientMagnitude.MagnitudeAutoCorrelation;
% Cross-Correlations
CrossMag0 = TexFeats.CoefficientMagnitude.OrientationCrossCorrelation;
CrossScalePhase0 = TexFeats.CoefficientMagnitude.ScaleCrossCorrelation;
% Cross-Correlation of Scales
CrossReal0 = TexFeats.CrossScalePhaseStatistics.ScaleCrossCorrelation;

% Check if the original image has the same dimensions as the requested
% Options.Size, in case ThisMask is not empty
ThisOriginal = imresize(ThisPyramid.Original,Options.Size); 
if nargin < 4, ThisMask = zeros(Options.Size); end;

% Initialize the Synthesis
ThisSynthesized = mu0 + sqrt(sigma0).*randn(Options.Size);
ThisPrevious = ThisSynthesized;

%% 2) Synthesizing Loop
ThisEvolution = zeros([size(ThisSynthesized,1),size(ThisSynthesized,2),size(ThisSynthesized,3),Options.MaxIter+1]);
ThisEvolution(:,:,:,1) = ThisSynthesized;

figure,
subplot(121), imshow(zeros(size(ThisSynthesized)),[]);
subplot(122), imshow(ThisSynthesized,[]);
for n = 1:Options.MaxIter
    
    fprintf('Synthesis Iteration: %i...\n',n);
    
    %% 2.1) Build Steerable Pyramid of the temp synthesized im
    fprintf('\tBuilding Steerable Pyramid...'); tic;
    SynthPyramid = buildSteerablePyramid(ThisSynthesized,ThisPyramid.Options,false);
    
    % subtract mean of lowband pyramid
    SynthPyramid.LowPass{end} = SynthPyramid.LowPass{end} - mean(SynthPyramid.LowPass{end}(:));
    
    % Extract Magnitude of the Pyramid bands.
    MagPyr = struct('HighPass',{abs(SynthPyramid.HighPass)},'LowPass',{cellfun(@(x) abs(x),SynthPyramid.LowPass,'un',0)},'Marginals',{cellfun(@(x) abs(x),SynthPyramid.Marginals,'un',0)});
    RePyr = struct('HighPass',{real(SynthPyramid.HighPass)},'LowPass',{cellfun(@(x) real(x),SynthPyramid.LowPass,'un',0)},'Marginals',{cellfun(@(x) real(x),SynthPyramid.Marginals,'un',0)});
    
    % Subtract mean from Magnitudes
    MuMag = struct('HighPass',{mean(MagPyr.HighPass(:))},'LowPass',{cellfun(@(x) cell2mat(arrayfun(@(n) mean(reshape(x(:,:,n),[],1)),permute((1:size(x,3)),[1 3 2]),'un',0)),MagPyr.LowPass,'un',0)},'Marginals',{cellfun(@(x) cell2mat(arrayfun(@(n) mean(reshape(x(:,:,n),[],1)),permute((1:size(x,3)),[1 3 2]),'un',0)),MagPyr.Marginals,'un',0)});
    MagPyr = struct('HighPass',{MagPyr.HighPass - MuMag.HighPass},'LowPass',{cellfun(@(x,y) (x-repmat(y,[size(x,1),size(x,2),1])),MagPyr.LowPass,MuMag.LowPass,'un',0)},'Marginals',{cellfun(@(x,y) (x-repmat(y,[size(x,1),size(x,2),1])),MagPyr.Marginals,MuMag.Marginals,'un',0)});
    
    fprintf('(%f sec)\n',toc);
    
    %% 2.2) Impose RawCoefficientCorrelation to synthesized LowBand of last
    % level
    fprintf('\tImposing Raw Coefficient Auto-Correlation to LowBand...'); tic;
    ThisSynthesized = RePyr.LowPass{end};
    options_tmp.NLevels = 0;
    options_tmp.NBands = 1;
    PyrTmp = buildSteerablePyramid(ThisSynthesized,options_tmp,false);
    ThisSynthesized = PyrTmp.LowPass{end};
    aCorr = rawCorr0{ThisPyramid.Options.NLevels+1};
    % If autocorr centered value is greater than 10^-4 times sigma0, then
    % impose autocorrelation; otherwise, simply add noise to the image:
    aCorr_center = aCorr(floor(WindowSize/2)+1,floor(WindowSize/2)+1);
    if (aCorr_center/sigma0) > 1e-4
        % Take only the part that we are interested in from aCorr (I'm not
        % really sure why Portilla et al did this on this part of the code).
        Sz = min(size(ThisSynthesized)/2);
        w1 = (WindowSize-1)/2;
        w0 = min(Sz/2-1,w1);
        aCorr = aCorr(w1-w0+1:w1+w0+1,w1-w0+1:w1+w0+1);
        ThisSynthesized = imposeAutoCorrelation(ThisSynthesized, aCorr);
    else
        ThisSynthesized = ThisSynthesized.*sqrt(aCorr_center/var(ThisSynthesized(:)));
    end
    % Discard (possible) imaginary part
    ThisSynthesized = real(ThisSynthesized);
    fprintf('(%f sec)\n',toc);
    
    %% 2.3) Impose Skewness 
    fprintf('\tImposing Skewness to LowBand...'); tic;
    ThisSynthesized = imposeSkewness(ThisSynthesized,skewLP(end));
    fprintf('(%f sec)\n',toc);
    
    %% 2.4) Impose Kurtosis
    fprintf('\tImposing Kurtosis to LowBand...'); tic;
    ThisSynthesized = imposeKurtosis(ThisSynthesized,kurtLP(end));
    fprintf('(%f sec)\n',toc);
        
    
    %% 2.5) Scale Reconstruction (in this loop we will reconstruct the 
    % result stored in ThisSynthesized (to which we previously imposed the
    % autocorrelation, the skewness and kurtosis), from coarser to fine
    % scale.
    for Scale = NLevels:-1:1
        
        fprintf('\tAnalyzing Level: %i\n',Scale);
        
        % Let's get the magnitudes of the pyramid at this scale (all
        % orientations), for the synthesized image. In the original code
        % this variable was called 'cousins'.
        HighLevelMags = MagPyr.LowPass{Scale}; 
        
        %% 2.5.1) Impose cross-Correlation of bands
        fprintf('\t\tImposing Cross-Correlation of Bands...'); tic;
        % To these magnitudes we want to impose the Cross-correlation
        % magnitude (OrientationCrossCorrelation) and, if possible, the 
        % ScaleCrossCorrelation of coarser levels (deeper levels). This
        % means that at the last scale (Scale = NLevels) we will impose
        % only OrientationCrossCorrelation, while at the rest of upper
        % scales we will impose both of them.
        if Scale ~= NLevels
            % impose both cross-correlation magnitude and
            % scalecrosscorrelaton
            % Let's get the lowlevels and expand them (and apply double
            % phase)
            LowLevelMags = SynthPyramid.LowPass{Scale+1};
            LowLevelMags = cell2mat(arrayfun(@(n) abs(dPhase(expand(LowLevelMags(:,:,n),2)/4)),permute((1:size(LowLevelMags,3)),[1 3 2]),'un',0));
            % subtract mean from lowlevelmags
            LowLevelMags = cell2mat(arrayfun(@(n) LowLevelMags(:,:,n) - sum(sum(LowLevelMags(:,:,n),1),2)/numel(LowLevelMags(:,:,n)),permute((1:size(LowLevelMags,3)),[1 3 2]),'un',0));
            HighLevelMags = imposeXCorr(HighLevelMags, CrossMag0{Scale}, LowLevelMags, CrossScalePhase0{Scale});
        else
            % impose only cross-correlation mag
            HighLevelMags = imposeXCorr(HighLevelMags, CrossMag0{Scale}, [], []);
        end
        fprintf('(%f sec)\n',toc);
        
        % Get only real part
        HighLevelMags = real(HighLevelMags);
        % Update Synthetic Pyramid
        MagPyr.LowPass{Scale} = HighLevelMags;

        
        %% 2.5.2) Impose auto Correlation of Magnitude Responses
        fprintf('\t\tImposing Auto-Correlation of Magnitude Responses...'); tic;
        aCorrMag = autoCrossMag0{Scale};
        MagTmp = MagPyr.LowPass{Scale};
        % Take only the part that we are interested in from aCorr (I'm not
        % really sure why Portilla et al did this on this part of the code).
        Sz = min([size(MagTmp,1),size(MagTmp,2)]/2);
        w1 = (WindowSize-1)/2;
        w0 = min(Sz/2-1,w1);
        aCorrMag = aCorrMag(w1-w0+1:w1+w0+1,w1-w0+1:w1+w0+1,:);
        % Now, for each orientation, imposeAutoCorrelation
        for b = 1:NBands
            subMagTmp = imposeAutoCorrelation(MagTmp(:,:,b), aCorrMag(:,:,b));
            subMagTmp = real(subMagTmp);
            % Put it back on the right place in MagPyr
            MagPyr.LowPass{Scale}(:,:,b) = subMagTmp;
            
            % Now Impose magnitude (discarding any negative number)
            subMagTmp = max(0,subMagTmp + magMeans0.LowPass{Scale}(:,:,b));
            % I don't understand the following line. I just copied it from 
            % the original code
            subMagTmp = subMagTmp./( abs(SynthPyramid.LowPass{Scale}(:,:,b)) + (abs(SynthPyramid.LowPass{Scale}(:,:,b)) < eps) );
            % Put it back on SynthPyramid
            SynthPyramid.LowPass{Scale}(:,:,b) = SynthPyramid.LowPass{Scale}(:,:,b).*subMagTmp;
        end
        fprintf('(%f sec)\n',toc);
        
        
        %% 2.5.3) Impose cross-correlation of Real Parts
        fprintf('\t\tImposing Cross-Correlation of Real Parts...'); tic;
        % Let's get the real part of the pyramid at this scale (all
        % orientations), for the synthesized image. In the original code
        % this variable was called 'cousins'.
        HighLevelReal = real(SynthPyramid.LowPass{Scale}); 
        % And the real part of the lowlevels (only if Scale ~= NLevels)
        if Scale ~= NLevels
            % Let's get the lowlevels and expand them (and apply double
            % phase)
            LowLevelReal = SynthPyramid.LowPass{Scale+1};
            LowLevelReal = cell2mat(arrayfun(@(n) dPhase(expand(LowLevelReal(:,:,n),2)/4),permute((1:size(LowLevelReal,3)),[1 3 2]),'un',0));
            % concatenate thru dimension 3
            LowLevelReal = cat(3,real(LowLevelReal),imag(LowLevelReal));
            
            % Now impose cross-correlation (real parts)
            HighLevelReal = cell2mat(arrayfun(@(n) imposeXCorr(HighLevelReal(:,:,n),mean(reshape(HighLevelReal(:,:,n).^2,[],1)),LowLevelReal,CrossReal0{Scale}(n,:)),permute((1:size(HighLevelReal,3)),[1 3 2]),'un',0));
            
        end
        % Update Synthetic pyramid
        SynthPyramid.LowPass{Scale} = HighLevelReal;
        fprintf('(%f sec)\n',toc);
        
        
        %% 2.5.4) Re-create analytic subbands
        fprintf('\t\tRecreating Analytic Subbands...'); tic;
        NewLowLevels = SynthPyramid.LowPass{Scale};
        % first create some auxiliar variables
        [M,N,K] = size(NewLowLevels);
        dims = [M,N];
        ctr = ceil((dims+0.5)/2);
        ang = mkAngle(dims, 0, ctr);
        ang(ctr(1),ctr(2)) = -pi/2;
        for b = 1:K
            % Apply phase to omega mask
            ang0 = pi*(b-1)/K;
            xang = mod(ang-ang0+pi, 2*pi) - pi;

            masktmp = 2*(abs(xang) < pi/2) + (abs(xang) == pi/2);
            masktmp(ctr(1),ctr(2)) = 1;
            masktmp(:,1) = 1;
            masktmp(1,:) = 1; 
            % Filter Data
            NewLowLevels(:,:,b) = ifft2(fftshift(masktmp).*fft2(NewLowLevels(:,:,b)));
        end
        % Update pyramid
        SynthPyramid.LowPass{Scale} = NewLowLevels;
        fprintf('(%f sec)\n',toc);
        
        
        %% 2.5.6) Reconstruction of Synthetic image
        fprintf('\t\tReconstucting Synthesized Image from Subbands...'); tic;
        % We have imposed all types of parameters to modify our noise into
        % what we want (texture). Now we must reconstruct the image using
        % all bands for this scale and subscales.
        
        % Now get only from this Scale to the end
        RePyrTmp = struct('HighPass',{zeros(size(SynthPyramid.LowPass{Scale}))},'LowPass',{cellfun(@(x) real(x),SynthPyramid.LowPass(Scale:Scale+1),'un',0)},'Marginals',{cellfun(@(x) real(x),SynthPyramid.Marginals(Scale:Scale+1),'un',0)});
        RePyrTmp.LowPass{end} = zeros([size(SynthPyramid.LowPass{Scale+1},1),size(SynthPyramid.LowPass{Scale+1},2)]);
        
        opts_tmp = ThisPyramid.Options;
        opts_tmp.NLevels = 1;
        ThisReconstruction = SteerableReconstruction(RePyrTmp, opts_tmp, false);
        
        % Now recombine previous Synthetic image with this Reconstruction
        ThisSynthesized = real(expand(ThisSynthesized,2))/4;
        ThisSynthesized = ThisSynthesized + ThisReconstruction{1};
        fprintf('(%f sec)\n',toc);
        
        
        %% 2.5.7) Impose Auto-correlation to reconstructed image
        fprintf('\t\tImposing Raw Coefficient Auto-Correlation to Reconstructed Image...'); tic;
        aCorr = rawCorr0{Scale};
        % If autocorr centered value is greater than 10^-4 times sigma0, then
        % impose autocorrelation; otherwise, simply add noise to the image:
        aCorr_center = aCorr(floor(WindowSize/2)+1,floor(WindowSize/2)+1);
        if (aCorr_center/sigma0) > 1e-4
            % Take only the part that we are interested in from aCorr (I'm not
            % really sure why Portilla et al did this on this part of the code).
            Sz = min(size(ThisSynthesized)/2);
            w1 = (WindowSize-1)/2;
            w0 = min(Sz/2-1,w1);
            aCorr = aCorr(w1-w0+1:w1+w0+1,w1-w0+1:w1+w0+1);
            ThisSynthesized = imposeAutoCorrelation(ThisSynthesized, aCorr);
        else
            ThisSynthesized = ThisSynthesized.*sqrt(aCorr_center/var(ThisSynthesized(:)));
        end
        ThisSynthesized = real(ThisSynthesized);
        fprintf('(%f sec)\n',toc);
        
        
        %% 2.5.8) Imposing Marginal Statistics to Reconstructed Image
        fprintf('\t\tImposing Skewness to Reconstructed Image...'); tic;
        ThisSynthesized = imposeSkewness(ThisSynthesized,skewLP(Scale));
        fprintf('(%f sec)\n',toc);
    
        fprintf('\t\tImposing Kurtosis to Reconstructed Image...'); tic;
        ThisSynthesized = imposeKurtosis(ThisSynthesized,kurtLP(Scale));
        fprintf('(%f sec)\n',toc);

    end
    
    
    %% 2.6) Impose HP Variance
    fprintf('\tImposing High-Frequency Variance to Reconstructed Image...'); tic;
    % We will only impose this condition in case the variance of the
    % reconstructed image (ThisSynthesized) is greater than the variance of
    % the original high-pass band (HPvar0)
    
    % Let's calculate, then, the variance of the high-frequency component
    % of the synthesized pyramid
    ThisHF = SynthPyramid.HighPass;
    HPvarX = var(ThisHF(:));
    if HPvarX > HPvar0
        ThisHF = ThisHF*sqrt(HPvar0/HPvarX);
        SynthPyramid.HighPass = ThisHF;

        [M,N,~] = size(ThisHF);

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
        rho = hypot(xgrid,ygrid);
        rho(centeri,centerj) = rho(centeri,centerj-1);
        rho = log2(rho);

        % 1.2) High pass mask (H0Mask) and high-pass component (H0). This component
        % will not be used in the building of the steerable pyramid.
        H0Mask = reshape(interp1(X,sqrt(Y),rho(:), 'linear', 'extrap'),size(rho));

        % apply fft 
        ThisHF = H0Mask.*fftshift(fft2(real(ThisHF)));
        % undo fft
        ThisHF = real(ifft2(ifftshift(ThisHF)));
    end
    
    % Add this reconstruction of HF to the synthesized image
    ThisSynthesized = ThisSynthesized + ThisHF;
    fprintf('(%f sec)\n',toc);

    %% 2.7) Impose Pixel Statistics
    fprintf('\tImposing Pixel Statistics to Reconstructed Image...'); tic;
    muX = mean(ThisSynthesized(:));
    varX = var(ThisSynthesized(:));
    % subtract mean
    ThisSynthesized = ThisSynthesized - muX;
    % adjust variance
    ThisSynthesized = ThisSynthesized.*sqrt(sigma0/varX);
    % adjust mean
    ThisSynthesized = ThisSynthesized + mu0;
    % impose skewness
    ThisSynthesized = imposeSkewness(ThisSynthesized,skew0);
    
    % impose kurtosis
    ThisSynthesized = imposeKurtosis(ThisSynthesized,kurt0);
    % impose range (simply cutting values above max0 or below min0). Notice
    % this it not like mat2gray!! We are not stretching the histogram of
    % ThisSynthesized, we're just cutting it!!
    ThisSynthesized = max(min(ThisSynthesized,max0),min0);
    %ThisSynthesized = ((max0-min0)/range(ThisSynthesized(:))).*(ThisSynthesized - min(ThisSynthesized(:))) + min0;
    fprintf('(%f sec)\n',toc);
    
    %% 2.9) Apply Mask. In case the user specified an input Mask we will 
    %  use it to put the original image part on ThisSynthesized at the 
    %  pixels where ThisMask == 1. This is used mainly in purposes of
    %  extrapolation of images (see Portilla's original paper, section 3.3
    fprintf('\tImposing Mask and Original image to Reconstructed Image...'); tic;
    ThisSynthesized = ThisMask.*ThisOriginal + (1-ThisMask).*ThisSynthesized;
    fprintf('(%f sec)\n',toc);
          
    %% 2.9) Apply Accelerator
    alpha = 0.8;
    Thistmp = ThisPrevious;
    ThisPrevious = ThisSynthesized;
    %% 2.9) Update figure
    D = ThisSynthesized - Thistmp;
    E0 = sum(abs(ThisSynthesized(:)));
    subplot(121), imshow(D,[]); title(sprintf('|D| = %3.2f%%',100*sum(abs(D(:)))/E0 ));
    subplot(122), imshow(ThisSynthesized,[]); title(sprintf('Iteration %i out of %i',n,Options.MaxIter));
    drawnow();
    
    ThisSynthesized = ThisSynthesized + alpha*(ThisSynthesized - Thistmp);
    
    ThisEvolution(:,:,:,n+1) = ThisPrevious;
    
    
    
end


ThisSynthesized = ThisPrevious;


end

function Options = checkInputOptions(TheseParameters)

OptionFields = {'MaxIter','Size'};
FieldsEmptyValues = {25;[256,256]};
Options = cell2struct(FieldsEmptyValues,OptionFields');

% Now fill data
if ~isempty(TheseParameters)
    if isstruct(TheseParameters)
        TheseParamsFields = fieldnames(TheseParameters);
        % Keep only fields that are valid
        TheseParamsFields = TheseParamsFields(cellfun(@(x) sum(strcmp(OptionFields,x)) > 0,TheseParamsFields));
        if ~isempty(TheseParamsFields)
           for i=1:numel(TheseParamsFields)
               Options.(TheseParamsFields{i}) = TheseParameters.(TheseParamsFields{i});
           end
        end
    end
end

% Make sure Size is a multiple of 2!!!
M = 2^floor(log2(Options.Size(1)));
N = 2^floor(log2(Options.Size(2)));

Options.Size = [M,N];

end

function Y = imposeAutoCorrelation(X, cY)

% Unnormalize previous autocorrelation cY
cY = cY.*numel(X);

% 1) Fourier Transform of X
[M,N,~] = size(X);
[K,~,~] = size(cY);
F = fft2(X);

% 2) Compute autocorrelation of X (cX)
autoCorr = @(X) fftshift(real(ifft2(abs(fft2(X)).^2)));
cX = autoCorr(X);
% Take the center part of cX
icenter = ceil((M-1)/2)+1;
jcenter = ceil((N-1)/2)+1;
w = ceil((K-1)/2); % width of window we are going to take
cX = cX(icenter-2*w:icenter+2*w,jcenter-2*w:jcenter+2*w);

% 3) Build the matrix that performs the convolution
wK = floor((K^2)/2) + 1;
TcX = im2col(cX,[K K],'sliding');
% Reshape and take only the part that we are interested in (1:wK)
TcX = reshape(arrayfun(@(n) reshape(TcX(:,n),[K K]),(1:size(TcX,2)),'un',0),[size(cX,1)-K+1 size(cX,2)-K+1])';
TcX = TcX(:);
TcX = TcX(1:wK);

factor = ones(K);
factor(wK) = 1/2;
TcX = cellfun(@(x) ((x + x(end:-1:1,end:-1:1)).*factor),TcX,'un',0);
TcX = cellfun(@(x) reshape(x',[],1),TcX,'un',0);
TcX = cellfun(@(x) (x(1:wK))',TcX,'un',0);
TcX = cell2mat(TcX);
% figure, imshow(cell2mat(TcX),[])

% 4) Take only the part required from cY (same as TcX)
cY = reshape(cY',[],1);
cY = cY(1:wK);

% 5) Calculate A_h (if you're lost look at Portilla et al. paper, addendum
% A, section A.2)
A_tmp = TcX\cY;
% And duplicate A_h (flipped)
A_tmp = cat(1,A_tmp,A_tmp(end-1:-1:1));
% And reshape to original size [K,K]
A_tmp = reshape(A_tmp,[K K])';

% Now expand A_h to [M,N] so we can convolve A_h and F
A_h = zeros(M,N);
A_h(icenter-w:icenter+w,jcenter-w:jcenter+w) = A_tmp;

% 6) Now take fourier transform of A_h
A_h = real(fft2(fftshift(A_h)));

% 7) Finally, convolve F with sqrt(abs(A_h)) to obtain new image with
% imposed autocorrelation
A_h = sqrt(abs(A_h));
Y = ifft2(F.*A_h);

end

function Y = imposeSkewness(X, skew0)

% 1) Subtract mean from input image
mu0 = mean(X(:));
X = X - mu0;


% 2) Extract moments from image (n=2 is variance of image, n=3 is related
% to skewness, etc.). For more info, take a look at this: 
% http://en.wikipedia.org/wiki/Moment_(mathematics)#Significance_of_the_moments
XMoments = permute(sum(sum(X.^repmat(permute((1:6),[1 3 2]),size(X)),1),2)./numel(X),[1 3 2]);
sdX = sqrt(XMoments(2));
skewX = XMoments(3)/sdX^3;


% 3) Compute coefficients of eta_n in equation (17) - See original paper.

% 3.1) Compute coefficients of the numerator p(1) + p(2)*lambda +
% p(3)*lambda^2 + p(4)*lambda^4
p = zeros(1,4);
p(4) = XMoments(6) - 3*sdX*skewX*XMoments(5) + 3*(sdX^2)*(skewX^2-1)*XMoments(4) + ...
       (sdX^6)*(2 + 3*skewX^2 - skewX^4);
p(3) = 3*(XMoments(5) - 2*sdX*skewX*XMoments(4) + (sdX^5)*(skewX^3));
p(2) = 3*(XMoments(4) - (sdX^4)*(1 + skewX^2));
p(1) = skewX*(sdX^3);

% 3.2) Compute coefficients of the denominator q(1) + q(2)*lambda +
% q(3)*lambda^3 (Please, notice that q(2) is actually zero, so won't be
% taken into account.
q = zeros(1,3);
q(3) = XMoments(4) - (1+skewX^2)*(sdX^4);
q(1) = sdX^2;

% 3.3) Find the coefficients 'a' from Equation (19).
a = zeros(1,7);
a(7) = p(4)^2;
a(6) = 2*p(4)*p(3);
a(5) = (p(3))^2 + 2*p(4)*p(2);
a(4) = 2*(p(4)*p(1) + p(3)*p(2));
a(3) = (p(2))^2 + 2*p(3)*p(1);
a(2) = 2*p(2)*p(1);
a(1) = (p(1))^2;

% 3.4) Find coefficients 'b' (It starts to get a lil bit messy here. See,
% in the original equation (19) you can see that some components have 2
% terms: the first term is calculated using coefficients from the numerator,
% like p(1),p(2), etc. while the second depends on the coefficients of the 
% denominator (q(1),q(3)). Instead of calculating simply a, we are dividing
% these two parts into 'a' and 'b'.
b = zeros(1,7);
b(7) = (q(3))^3;
b(5) = 3*q(1)*(q(3)^2);
b(3) = 3*q(3)*(q(1))^2;
b(1) = (q(1))^3;


% 4) Now Compute derivative with respect to lambda 
d = zeros(1,8);
d(8) = p(3)*b(7);
d(7) = 2*p(2)*b(7) - p(4)*b(5);
d(6) = 3*p(1)*b(7);
d(5) = p(2)*b(5) - 2*p(4)*b(3);
d(4) = 2*p(1)*b(5) - p(3)*b(3);
d(3) = -3*p(4)*b(1);
d(2) = p(1)*b(3) - 2*p(3)*b(1);
d(1) = -p(2)*b(1);

% calculate roots of derivative (solutions of lambda)
Lambda = roots(d(end:-1:1));

% Now let's keep only the solutions that do not (or almost don't) have
% imaginary part. In order to do this we calculate first the phase of the
% solutions (imag(Lambda)/real(Lambda)) and check if that number of close
% to zero
Lambda = real(Lambda(abs(imag(Lambda)./real(Lambda)) < 1e6));

% We need to find the largest and smallest Lambda so the maximum positive
% Root and the maximum negative root too.
Lambda = [-min(abs(min(Lambda(Lambda<0))),1/eps),... 
           min(max(Lambda(Lambda>=0)),1/eps)];


% 5) Using these roots let's calculate the minimum and maximum skewness to
% be imposed to the image. 
skewY = polyval(p(end:-1:1),Lambda)./(polyval(b(end:-1:1),Lambda)).^0.5;
% The skewness to be imposed (skew0) should be between skewY range
if skew0 <= min(skewY) || skew0 >= max(skewY)
    Lambda = max(min(max(Lambda),skew0),min(Lambda));
    
    % adjust skewness
    Y = X + Lambda*(X.^2 - sdX.^2 -sdX*skewX.*X );
    % adjust variance
    Y = Y*sqrt(XMoments(2)/mean(Y(:).^2));
    %adjust mean
    Y = Y + mu0;
    
else
    
    % Find coefficients c
    c = a - b.*(skew0^2);
    % find roots of c (solutions to the poly)
    r = roots(c(end:-1:1));
    
    % Keep only roots that have almost no imaginary part (like we did
    % before)
    r = r(abs(imag(r)./real(r)) < 1e-6);
    % Keep only roots whose real part has the same sign as (skew0 - skewX)
    r = real(r(sign(r) == sign(skew0-skewX)));

    % if no root matched the previous conditions, return and skip skewness
    % adjustment
    if isempty(r)
        Y = X*sqrt(XMoments(2)/mean(X(:).^2));
        %adjust mean
        Y = Y + mu0;
        warning('\tSkipping Skewness adjustment...\n');
        return;
    end
    
    % If we have made it so far, let's find the values of the r (solve the
    % system)
    rsign = sign(polyval(p(end:-1:1),r));
    % Let's discard any value of r that has different sign than skew0 or 0
    r = r((rsign == 0)  | (rsign == sign(skew0)));
    
    % if no root matched the previous conditions, return and skip skewness
    % adjustment
    if isempty(r)
        Y = X*sqrt(XMoments(2)/mean(X(:).^2));
        %adjust mean
        Y = Y + mu0;
        warning('\tSkipping Skewness adjustment...\n');
        return;
    end
    
    % If r is not empty let's get the smallest value of r that fix the skew
    Lambda = r(find(abs(r) == min(abs(r)),1));
    
    % adjust skewness
    Y = X + Lambda*(X.^2 - sdX.^2 -sdX*skewX.*X );
    % adjust variance
    Y = Y*sqrt(XMoments(2)/mean(Y(:).^2));
    %adjust mean
    Y = Y + mu0;

end


end

function Y = imposeKurtosis(X, kurt0)

% 1) Subtract mean from input image
mu0 = mean(X(:));
X = X - mu0;


% 2) Extract moments from image (n=2 is variance of image, n=3 is related
% to skewness, etc.). For more info, take a look at this: 
% http://en.wikipedia.org/wiki/Moment_(mathematics)#Significance_of_the_moments
XMoments = permute(sum(sum(X.^repmat(permute((1:12),[1 3 2]),size(X)),1),2)./numel(X),[1 3 2]);
%kurtX = XMoments(4)/XMoments(2)^2;
alpha = XMoments(4)/XMoments(2);


% 3) Find the coefficients of equation (21)
% 3.1) Numerator
p = zeros(1,5);
p(5) = XMoments(12)-4*alpha*XMoments(10)-4*XMoments(3)*XMoments(9)+6*alpha^2*XMoments(8)+12*alpha*XMoments(3)*XMoments(7)+6*XMoments(3)^2*XMoments(6)-...
	   4*alpha^3*XMoments(6)-12*alpha^2*XMoments(3)*XMoments(5)+alpha^4*XMoments(4)-12*alpha*XMoments(3)^2*XMoments(4)+...
	   4*alpha^3*XMoments(3)^2+6*alpha^2*XMoments(3)^2*XMoments(2)-3*XMoments(3)^4;
p(4) = 4*(XMoments(10)-3*alpha*XMoments(8)-3*XMoments(3)*XMoments(7)+3*alpha^2*XMoments(6)+6*alpha*XMoments(3)*XMoments(5)+3*XMoments(3)^2*XMoments(4)-...
	   alpha^3*XMoments(4)-3*alpha^2*XMoments(3)^2-3*XMoments(4)*XMoments(3)^2);
p(3) = 6*(XMoments(8)-2*alpha*XMoments(6)-2*XMoments(3)*XMoments(5)+alpha^2*XMoments(4)+2*alpha*XMoments(3)^2+XMoments(3)^2*XMoments(2));
p(2) = 4*(XMoments(6)-alpha^2*XMoments(2)-XMoments(3)^2);
p(1) = XMoments(4);

% 3.2) Denominator
q = zeros(1,3);
q(1) = XMoments(2);
q(3) = p(2)/4;


% 4) Compute derivative with respect to lambda
d = zeros(1,5);
d(1) = p(4)*q(3);
d(2) = 2*p(3)*q(3) - 4*p(5)*q(1);
d(3) = 4*q(3)*p(2) -3*p(4)*q(1) - p(2)*q(3);
d(4) = 4*q(3)*p(1) - 2*p(3)*q(1);
d(5) = -p(2)*q(1);


% calculate roots of derivative (solutions of lambda)
Lambda = roots(d);

% Now let's keep only the solutions that do not (or almost don't) have
% imaginary part. In order to do this we calculate first the phase of the
% solutions (imag(Lambda)/real(Lambda)) and check if that number of close
% to zero
Lambda = real(Lambda(abs(imag(Lambda)./real(Lambda)) < 1e-6));

LNeg = Lambda(Lambda < 0);
LPos = Lambda(Lambda >= 0);

% We need to find the largest and smallest Lambda so the maximum positive
% Root and the maximum negative root too.
Lambda = [-min(abs(min(Lambda(Lambda<0))),1/eps),... 
           min(max(Lambda(Lambda>=0)),1/eps)];

% 5) Find solutions to equation (21) for Lambdas
kurtY = polyval(p(end:-1:1),Lambda)./polyval(q(end:-1:1),Lambda).^2;
% The kurtosis to be imposed (kurt0) should be between kurtY range
if kurt0 <= min(kurtY) 
    Lambda = max(LNeg);
    
    Y = X + Lambda*(X.^3 - alpha.*X - XMoments(3) );
    % adjust variance
    Y = Y*sqrt(XMoments(2)/var(Y(:)));
    %adjust mean
    Y = Y + mu0;
    
elseif kurt0 >= max(kurtY)
    Lambda = min(LPos);
   
    % adjust kurtosis
    Y = X + Lambda*(X.^3 - alpha.*X - XMoments(3) );
    % adjust variance
    Y = Y*sqrt(XMoments(2)/var(Y(:)));
    %adjust mean
    Y = Y + mu0;
    
else
    
    % find the coefficients c to solve the system
    c = zeros(1,5);
    c(5) = p(1) - kurt0*q(1)^2;
    c(4) = p(2);
    c(3) = p(3) - 2*kurt0*q(3)*q(1);
    c(2) = p(4);
    c(1) = p(5) - kurt0*q(3)^2;
    
    % get the roots
    r = roots(c);
    
    % Keep only roots that have almost no imaginary part (like we did
    % before)
    r = real(r(abs(imag(r)./real(r)) < 1e-6));
    
    % if no root matched the previous conditions, return and skip skewness
    % adjustment
    if isempty(r)
        Y = X*sqrt(XMoments(2)/var(Y(:)));
        %adjust mean
        Y = Y + mu0;
        warning('\tSkipping Kurtosis adjustment...\n');
        return;
    end
    
    % If r is not empty let's get the smallest value of r that fix the skew
    Lambda = r(find(abs(r) == min(abs(r)),1));
    
    % adjust skewness
    Y = X + Lambda*(X.^3 - alpha.*X - XMoments(3) );
    % adjust variance
    Y = Y*sqrt(XMoments(2)/var(Y(:)));
    %adjust mean
    Y = Y + mu0;
    
end


end

function Z = imposeXCorr(HighLevel0, CrossMag0, LowLevel0, CrossScale0)

% Depending on whether we have both HighLevel0 and LowLevel0 we will impose
% only CrossMag0 or both CrossMag0 and CrossScale0. Let's start with
% CrossMag0.

% Check which inputs we have
if isempty(HighLevel0) || isempty(CrossMag0)
    Z = HighLevel0;
    fprintf('\tSkipping Cross-Correlation adjustment...\n');
    return;
end
    
[M,N,K] = size(HighLevel0);

% First let's vectorize HighLevel0 to carry on the arithmetic
% operations we must apply
X = reshape(HighLevel0,M*N,K,1);

if isempty(LowLevel0) || isempty(CrossScale0)  % Apply only CrossMag0

    % Now let's calculate the cross-correlation of the synthesized image (C
    % is the cross-correlation between bands) -> Equation (34) and previous
    % paragraph
    C = X'*X / (M*N);
    % Now compute eigenvectors of cross-correlation matrix (Equation (35))
    [E, D] = eig(C);
    D = diag(D);
    % Sort from greater to lower
    [D,idx] = sort(D,'descend');
    E = E(:,idx);
    % Reconstruct the diagonal matrix using D
    D = diag(sqrt(D));

    % 2) Now let's apply the same methodology (eigenvectors...) to input
    % cross-correlation
    [E0,D0] = eig(CrossMag0);
    D0 = diag(D0);
    [D0,idx] = sort(D0,'descend');
    E0 = E0(:,idx);
    D0 = diag(sqrt(D0));

    % 3) Now let's find the orthogonal matrix O: Equation (37). This Matrix
    % relates the cross-correlations of the synthesized image (E) and the
    % desired cross-correlation we want to impose (E0)
    O = E'*E0;

    % 4) Finally, the matrix M that we must apply to the input HighLevels0
    % in order to impose these CrossMag0 is: Equation (36)
    M = (E/D)*O*D0*E0';

    % 5) Apply M matrix transformation to input images (vectorized)
    X = X*M;
    % And reshape back
    Z = reshape(X,size(HighLevel0));

else % Apply both

    % Now let's calculate the cross-correlation of the synthesized image (C
    % is the cross-correlation between bands) -> Equation (34) and previous
    % paragraph
    Bx = X'*X / (M*N);
    
    [M,N,K] = size(LowLevel0);

    % First let's vectorize HighLevel0 to carry on the arithmetic
    % operations we must apply
    Y = reshape(LowLevel0,M*N,K,1);
    Bxy = (X'*Y)/(M*N);
    By = Y'*Y / (M*N);
    
    % The Synthesized M is:
    Msyn = Bx - ((Bxy/By) * Bxy');
    % While the input M is:
    Min = CrossMag0 - ((CrossScale0/By)*CrossScale0');
    
    % Now compute eigenvectors of cross-correlation matrix (Equation (35))
    [E, D] = eig(Msyn);
    D = diag(D);
    % Sort from greater to lower
    [D,idx] = sort(D,'descend');
    E = E(:,idx);
    % Reconstruct the diagonal matrix using D
    D = diag(sqrt(D));

    % 2) Now let's apply the same methodology (eigenvectors...) to input
    % cross-correlation
    [E0,D0] = eig(Min);
    D0 = diag(D0);
    [D0,idx] = sort(D0,'descend');
    E0 = E0(:,idx);
    D0 = diag(sqrt(D0));

    % 3) Now let's find the orthogonal matrix O: Equation (37). This Matrix
    % relates the cross-correlations of the synthesized image (E) and the
    % desired cross-correlation we want to impose (E0)
    O = E'*E0;
    
    % 4) Find the matrices Mx and My of transformation (for each type of
    % cross-correlation)
    Mx =  (E/D)*O*D0*E0';
    My =  By\(CrossScale0' - Bxy' * Mx);
    
    % 5) Apply transformation Mx and My to input values (this is the
    % imposition and adjustment)
    Z = X*Mx + Y*My;
    % And reshape back
    Z = reshape(Z,size(HighLevel0));
end


end