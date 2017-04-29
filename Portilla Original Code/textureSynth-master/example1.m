% Example 1: Synthesis of a "text" texture image, using
% Portilla-Simoncelli texture analysis/synthesis code, based on
% alternate projections onto statistical constraints in a complex
% overcomplete wavelet representation.
%
% See Readme.txt, and headers of textureAnalysis.m and
% textureSynthesis.m for more details.
%
% Javier Portilla (javier@decsai.ugr.es).  March, 2001

close all

im0 = pgmRead('text.pgm');	% im0 is a double float matrix!

files = dir('E:\DATASETS\KTHTIPS\*.png');
files = {files.name};

im0 = imresize(double(imread(fullfile('E:\DATASETS\KTHTIPS',files{randi(numel(files))}))),size(im0));
% k = randi(size(tt,1));
% im0 = imresize(tt(k:k+1000,:),size(im0));
% im0 = tt;
% sz = 2.^floor(log2(size(im0)));
% im0 = tt(1:sz(1),1:sz(2));
Nsc = 4; % Number of scales
Nor = 4; % Number of orientations
Na = 9;  % Spatial neighborhood is Na x Na coefficients
	 % It must be an odd number!

% params = texture_synthesis(im0);
     
params = textureAnalysis(im0, Nsc, Nor, Na);

Niter = 25;	% Number of iterations of synthesis loop
Nsx = size(im0,1);	% Size of synthetic image is Nsy x Nsx
Nsy = size(im0,2);	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

res = textureSynthesis(params, [Nsy Nsx], Niter);


figure;
subplot(121)
showIm(im0, 'auto', 1, 'Original texture');
subplot(122);
showIm(res, 'auto', 1, 'Synthesized texture');

% Can you read the NEW text? ;-)
