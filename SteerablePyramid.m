function [ThisPyramid,TexFeats] = SteerablePyramid(ThisImage,Options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of Steerable Pyramid for texture characterization & Synthesis
% as explained by Javier Portilla and Eero Simoncelli.
% Work described in:
%  "A Parametric Texture Model based on Joint Statistics of Complex Wavelet 
%   Coefficients".
%  J Portilla and E P Simoncelli. Int'l Journal of Computer Vision,
%  vol.40(1), pp. 49-71, Dec 2000.
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

if nargin == 0
    error('This function requires, at least, one image as input.');
elseif nargin == 1
    Options = [];
    Options.NLevels = 4; % Number of levels of the pyramid
    Options.NBands = 4; % Number of orientations for each level of the pyramid
    Options.WindowSize = 9; % Window size of the processing block (9 px by 9 px)
end

if ~ismatrix(ThisImage) || isempty(ThisImage)
    error('This function only accepts non-empty gray matrix images as inputs.');
end

% Check options struct
[ThisImage,Options] = CheckInputs(ThisImage,Options);

%% Add noise to the image
noise = range(ThisImage(:))/1000*randn(size(ThisImage));
ThisImage = ThisImage + noise;

%% 1) Building of the Steerable Pyramid
ThisPyramid = buildSteerablePyramid(ThisImage,Options);

% Add parameters to ThisPyramid before leaving
ThisPyramid.Options = Options;
ThisPyramid.Original = ThisImage - noise;

%% 2) Extraction of textural parameters that will characterize the image 
% (using the pyramid) much better than just TexFeats
% do this only if the user asked for it
if nargout == 1, TexFeats = []; return; end;
TexFeats = SteerableTexturalParams(ThisImage, ThisPyramid,Options);



end


% Checking the Inputs function
function [ThisImage,Options] = CheckInputs(ThisImage,TheseParameters)

OptionFields = {'NLevels','NBands','WindowSize'};
FieldsEmptyValues = {4;4;9};
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

% Now Check some values: 

% let's check the size of the image. We have to round its size to make
% sure it has log2 size (for instance, if we have an image whose size is
% [310, 1420] we should rescale it to [256,1024], otherwise when we
% construct our steerable pyramid we won't be able to successfully divide
% our image by 2 ([310, 1420]/2 -> [155, 710]/2 -> [77.5, 355]/2 -> etc)
[M,N,~] = size(ThisImage);
Options.OriginalImageSize = [M,N];

M = 2^floor(log2(M));
N = 2^floor(log2(N));

ThisImage = imresize(ThisImage,[M N]);

% 1) WindowSize must be an odd number, greater than 1
if Options.WindowSize <= 1
    Options.WindowSize = 3;
    warning('Window Size must be an odd number greater than 1. Setting its value to 3 and continuing.');
end
% round to closest odd number
Options.WindowSize = 2*round((Options.WindowSize+1)/2)-1;

% 2) Number of levels must be equal or greater than 1 and lower than
% floor(log2(min(M,N))) - 2
if isempty(Options.NLevels), Options.NLevels = floor(log2(min(M,N)/Options.WindowSize)); end;

if Options.NLevels < 1
    Options.NLevels = 1;
    warning('Number of Levels must be equal or greater than 1. Setting value to 1 and continuing.');
elseif Options.NLevels > floor(log2(min(M,N)/Options.WindowSize))
    Options.NLevels = floor(log2(min(M,N)/Options.WindowSize));
    warning('Number of Levels must be equal or smaller than %i. Setting value to %i and continuing.',Options.NLevels,Options.NLevels);
end


end


% Extraction of extra parameters for complete texture Characterization
function TexFeats = SteerableTexturalParams(ThisImage, ThisPyramid, Options)

fprintf('Extracting Features for Texture Characterization...\n'); 
% We will use both the magnitude of the pyramid and the real part of it:
MagPyr = struct('HighPass',{abs(ThisPyramid.HighPass)},'LowPass',{cellfun(@(x) abs(x),ThisPyramid.LowPass,'un',0)},'Marginals',{cellfun(@(x) abs(x),ThisPyramid.Marginals,'un',0)});
RePyr = struct('HighPass',{real(ThisPyramid.HighPass)},'LowPass',{cellfun(@(x) real(x),ThisPyramid.LowPass,'un',0)},'Marginals',{cellfun(@(x) real(x),ThisPyramid.Marginals,'un',0)});

% First thing, subtract mean of magnitude measures (and store this value):
muMag = struct('HighPass',{mean(mean(MagPyr.HighPass,1),2)},'LowPass',{cellfun(@(x) mean(mean(x,1),2),MagPyr.LowPass,'un',0)},'Marginals',{cellfun(@(x) mean(mean(x,1),2),MagPyr.Marginals,'un',0)});
MagPyr = struct('HighPass',{MagPyr.HighPass - muMag.HighPass},'LowPass',{cellfun(@(x,y) (x-y),MagPyr.LowPass,muMag.LowPass,'un',0)},'Marginals',{cellfun(@(x,y) (x-y),MagPyr.Marginals,muMag.Marginals,'un',0)});

% Now, let's reconstruct the images using the steerable pyramid (kinda like
% the opposite of constructing the pyramid)
% ThisReconstruction = SteerableReconstruction(RePyr, Options, false);

% Build Pyramid using lowband
opts_tmp.NLevels = 0;
opts_tmp.NBands = 1;
PyrTmp = buildSteerablePyramid(real(ThisPyramid.LowPass{end}),opts_tmp,false);

ThisReconstruction{Options.NLevels+1} = PyrTmp.LowPass{1};
for n=Options.NLevels:-1:1
    subPyrTmp = RePyr.LowPass{n};
    tnull = zeros(size(subPyrTmp,1),size(subPyrTmp,2));
    tnull2 = imresize(zeros(size(subPyrTmp,1),size(subPyrTmp,2)),0.5);
    subPyrTmp = struct('HighPass',{tnull},...
                       'LowPass',{[{subPyrTmp},{tnull2}]},...
                       'Marginals',{tnull});
    opts_tmp.NLevels = 1;
    opts_tmp.NBands = 4;               
    imTmp = SteerableReconstruction(subPyrTmp,opts_tmp,false);
    imTmp = imTmp{1};
    
    imReconstructed = real(expand(ThisReconstruction{n+1},2))/4;
    ThisReconstruction{n} = imReconstructed + imTmp;
    
end

% Every calculation on lowband should be done using PyrTmp.LowPass{1}

fprintf('\tMarginal Statistics...'); tic;
%% a) Marginal Statistics (TexFeats.MarginalStatistics)
% Opt0.NLevels = 0;
% Opt0.NBands = 0;
% Opt0.OriginalImageSize = size(ThisPyramid.LowPass{end});
% Pyr0 = buildSteerablePyramid(ThisPyramid.LowPass{end},Opt0);

% a.1) skewness and kurtosis of the partially reconstructed lowpass images 
% at each scale (2*(NLevels + 1) parameters) -> NLevels +1 because we are
% going to compute the skewness of the reconstructed images. In order to
% reconstruct the images what we are going to do is the following: Start by
% the lowest level (ThisPyramid.LowPass{end}), combine all lowpass images,
% and upsample, then repeat this procedure until you reach the top level.
% That way you have undone the process that constructed the Steering
% Pyramid. Each one of these Levels will produce a partial reconstruction.
% From each one of these partial reconstructions we can compute all
% parameters we need (autocorrelation, skewness, etc.). The 'plus 1' in
% (NLevels + 1) is due to the fact that our pyramid has a marginal image,
% which is what's left of the image after all downscaling and filtering
% processes.
% In the original code this feature was called statsLPim.
% kurtosis
kurtRecon = cellfun(@(x) kurtosis(x(:)),ThisReconstruction,'un',1);
% skewness
skewRecon = cellfun(@(x) skewness(x(:)),ThisReconstruction,'un',1);
TexFeats.MarginalStatistics.LP = struct('Kurtosis',kurtRecon,'Skewness',skewRecon);

% a.2) variance of the high-pass band (H0)(1 parameter). In the original
% code this feature was called: 'varianceHPR'
TexFeats.MarginalStatistics.HP = var(ThisPyramid.HighPass(:));

% a.3) mean, variance, skew, kurtosis, minimum and maximum values of the 
%      image pixels (original image) (6 parameters). In the original code
%      this feature was called 'statg0'
TexFeats.MarginalStatistics.Original = [mean(ThisImage(:)), var(ThisImage(:)), skewness(ThisImage(:)),...
            kurtosis(ThisImage(:)), min(ThisImage(:)), max(ThisImage(:))];

% a.4) mean of magnitudes. In the original code this feature was called
% 'magMeans0'
TexFeats.MarginalStatistics.MagnitudeMeans = muMag;

fprintf('(%f sec)\n',toc);

fprintf('\tRaw Coefficient Correlation...');
%% b) Raw Coefficient Correlation (TexFeats.RawCoefficientCorrelation)
% autocorrelation is defined as the convolution of a certain function (f)
% with itself. In fourier domain this can be translated as a product of the
% transform of f -> F(w), squared: autoCorrelation(f) = fourier(f) conv
% fourier(f) = (fourier(f))^2. 
autoCorr = @(X) fftshift(real(ifft2(abs(fft2(X)).^2)))/numel(X);
centeredSample = @(X,WindowSize)   X( (size(X,1)/2 + 1 -(WindowSize-1)/2):(size(X,1)/2 + 1 +(WindowSize-1)/2), (size(X,2)/2 + 1 -(WindowSize-1)/2):(size(X,2)/2 + 1 +(WindowSize-1)/2), :, :);
% auto correlation of real part of partially reconstructed images. In the
% original code this feature was called 'acr'
autoCorrReal = cellfun(@(x) centeredSample(autoCorr(x),Options.WindowSize),ThisReconstruction,'un',0);
% add lowpass remainings
autoCorrReal{end+1} = centeredSample(autoCorr(PyrTmp.LowPass{1}),Options.WindowSize);
TexFeats.RawCoefficientCorrelation = autoCorrReal;

fprintf('(%f sec)\n',toc);

fprintf('\tCoefficient Magnitude...');
%% c) Coefficient Magnitude (TexFeats.CoefficientMagnitude)
% c.1. MagnitudeAutoCorrelation. In the original code this feature was called: 'ace'
autoCorrMag = cellfun(@(x) centeredSample(cell2mat(arrayfun(@(n) autoCorr(x(:,:,n)),permute((1:size(x,3)),[1 3 2]),'un',0)),Options.WindowSize),MagPyr.LowPass(1:end-1),'un',0);
TexFeats.CoefficientMagnitude.MagnitudeAutoCorrelation = autoCorrMag;

% In order to extract all following parameters we need to compute what the
% authors of the paper (portilla et al.) call 'double-phased' components.
% The formula to obtain this component by using the real and imaginary
% parts of a matrix M is:
% doublePhased = abs(M).*exp( 2*(1i)*atan2(real(M),imag(M) )
dPhase = @(M) abs(M).*exp( (2i).*atan2(real(M),imag(M) ) );
%doublePhased = struct('HighPass',{dPhase(ThisPyramid.HighPass)},'LowPass',{cellfun(@(x) dPhase(x),ThisPyramid.LowPass,'un',0)},'Marginals',{cellfun(@(x) dPhase(x),ThisPyramid.Marginals,'un',0)});


% c.2. OrientationCrossCorrelation (cross-correlation of each subband
% magnitudes with those of other orientations at the same scale). In the
% original code this feature was called 'cousinMagCorr'.
A = cellfun(@(x) reshape(x,[],size(x,3),1),MagPyr.LowPass,'un',0);
OrientationCrossCorrelation = cellfun(@(x) (x'*x)./size(x,1),A,'un',0);
TexFeats.CoefficientMagnitude.OrientationCrossCorrelation = OrientationCrossCorrelation;

% c.3. ScaleCrossCorrelation (cross-correlation of subband magnitudes
% with all orientations at a coarser scale). In the original code this
% feature was called 'parentMagCorr'.
%B = cellfun(@(x) cell2mat(arrayfun(@(n) expand(x(:,:,n),2)/4,permute((1:size(x,3)),[1 3 2]),'un',0)),doublePhased.LowPass(2:end-1),'un',0);
B = cellfun(@(x) dPhase(cell2mat(arrayfun(@(n) expand(x(:,:,n),2)/4,permute((1:size(x,3)),[1 3 2]),'un',0))),ThisPyramid.LowPass(2:end-1),'un',0);
B = cellfun(@(x) reshape(abs(x),[],size(x,3),1),B,'un',0);
% subtract mean
B = cellfun(@(x) bsxfun(@minus,x,mean(x,1)),B,'un',0);
ScaleCrossCorrelation = cellfun(@(x,y) (x'*y)./size(x,1),A(1:end-2),B,'un',0);
TexFeats.CoefficientMagnitude.ScaleCrossCorrelation = ScaleCrossCorrelation;

fprintf('(%f sec)\n',toc);

fprintf('\tCross-scale phase statistics...');
%% d) Cross-scale phase statistics (TexFeats.CrossScalePhaseStatistics)
% d.1. OrientationCrossPhase (like OrientationCrossCorrelation, but with
% real parts instead of magnitudes). In the original code this feature was
% called 'cousinRealCorr'.
A = cellfun(@(x) reshape(x,[],size(x,3),1),RePyr.LowPass,'un',0);
OrientationCrossPhase = cellfun(@(x) (x'*x)./size(x,1),A,'un',0);
TexFeats.CrossScalePhaseStatistics.OrientationCrossCorrelation = OrientationCrossPhase;

% d.2. ScaleCrossCorrelation (cross-correlation of the real part of
% coefficients with both the real and imaginary part of the phase-doubled
% coefficients at all orientations at the next coarser scale). In the
% original code this feature was called: 'parentRealCorr'.
B = cellfun(@(x) dPhase(cell2mat(arrayfun(@(n) expand(x(:,:,n),2)/4,permute((1:size(x,3)),[1 3 2]),'un',0))),ThisPyramid.LowPass(2:end-1),'un',0);
B = cellfun(@(x) cat(2,reshape(real(x),[],size(x,3),1),reshape(imag(x),[],size(x,3),1)),B,'un',0);
ScaleCrossPhase = cellfun(@(x,y) (x'*y)./size(x,1),A(1:end-2),B,'un',0);
TexFeats.CrossScalePhaseStatistics.ScaleCrossCorrelation = ScaleCrossPhase;

fprintf('(%f sec)\n',toc);


end



