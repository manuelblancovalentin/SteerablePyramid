%% DEMONSTRATION OF PORTILLA'S STEERABLE PYRAMID TECHNIQUE FOR TEXTURE
%  CHARACTERIZATION AND SYNTHESIS. 
% 
% This code is based on Portilla's and Simoncelli's paper:
%  "A Parametric Texture Model based on Joint Statistics of Complex Wavelet 
%   Coefficients".
%  J Portilla and E P Simoncelli. Int'l Journal of Computer Vision,
%  vol.40(1), pp. 49-71, Dec 2000.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Manuel Blanco Valentin
% Industrial Engineer & Project Analyst 
% (Polytechnical University of Catalonia - EUETIB - Barcelona - Spain)
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
% 
%  I created this code for my MSc. Thesis in Oil Reservoir Image Data
%  processing as I needed to understand how the workflow proposed by
%  Portilla et al worked. Even though Portilla's original code (see here:
%  http://www.cns.nyu.edu/~eero/steerpyr/) works great, it was very
%  difficult to fully comprehend each step in the different processes of
%  building the steerable pyramid, characterizing the textures and
%  synthesizing them, as almost no comments nor intelligible variable names 
%  were used on the code. 
% 
%  I studied the code for several weeks and implemented all processes on my
%  own. I commented all steps and tried to include as many references to
%  the original paper as possible; so If anyone wants/needs to
%  understand how this code works, it'll be a lil' bit easier. 
% 
%  I have implemented ALL functions from scratch and I only use 1 of
%  Portilla's original code (expand function). I have also included another
%  function (dispPyramid) which displays the steerable pyramid in a more
%  friendly way (giving the impression of a real pyramid), so that if you
%  need to display the pyramid in any paper or work, you can simply use
%  this function almost out the box.
% 
%  I have tested my code and Portilla's original code using the same input
%  images and same initial conditions (starting white noise) and they
%  provide the EXACT same results (I checked this value per value).
% 
%  This code can be used mainly for three purposes: 
% 
%   1) Build a steerable pyramid, given an input image, and a number of
%   desired scales and orientations. A steerable pyramid is a technique
%   that uses a recursive filtering workflow to obtain the different
%   attributes (texturally speaking) that an image may have at different
%   scales and orientations. This workflow basically consists on taking the
%   input image, filtering it with a highpass filter, a series of
%   different-oriented 2D band-pass filters and a low-pass filter. The
%   low-pass component (what's left-over) is then downscaled using a factor
%   of 2 and all the process is repeated again (High-pass + band-pass +
%   LP). By doing so we are extracting textural information about
%   orientation of the features in our image at different scales. 
% 
%   2) Extract Features that characterize the texture on the input image at
%   different scales and orientations (based on the steerable pyramid built
%   previously). These feature are those described on Portilla's paper too.
%   As shown on that paper they seem to be sufficient to characterize
%   several types of textures well (indeed they seem to do its job so well
%   that we can actually synthesize almost identical textures by using
%   them, as I will explain on the following entry). On the other side,
%   Portilla and Simoncelli asure they are universal parameters (which
%   roughly means that you can use it in a bunch of different applications
%   without having to implement much changes).
% 
%   3) Texture Synthesis. As described in Portilla's paper, we can use the
%   textural features extracted from the steerable pyramid (used to
%   characterize it) to create a totally artificial and synthetic image
%   that will have the same textural characteristics as the original one.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HOW TO USE THIS CODE

% 1) Read the image you want to characterize
addpath(genpath(fullfile(pwd,'test-images')));
files = dir([fullfile(pwd,'test-images'),'\*.pgm']);
files = {files.name};

im0 = pgmRead(files{randi(numel(files))});
figure, imshow(im0,[])

% 2) Build the Steerable Pyramid and extract the textural features from it,
% specifiying the number of levels and bands, as well as the window size of
% the sliding window used for texture feature characterization. 
options.NLevels = []; % Number of levels of the pyramid. Leave empty for maximum allowed number of levels
options.NBands = 4; % Number of orientations for each level of the pyramid
options.WindowSize = 9; % Window size of the processing block (9 px by 9 px)

% Build pyramid
[SPyr,TexFeats] = SteerablePyramid(im0,options);

% Let's take a look at the Steerable Pyramid we just built
dispPyramid(im0,SPyr.LowPass(1:end-1),SPyr.Marginals(1:end-1))

% Let's say you want to take a look at a certain band (3) and scale (1) of 
% the pyramid. The syntaxis is like is:
figure,
imshow(SPyr.LowPass{1}(:,:,3),[])

% We can also reconstruct the original image using the pyramid we just
% created
ThisReconstruction = SteerableReconstruction(SPyr, SPyr.Options);
% Let's compare it with the original image
figure, 
subplot(121), imshow(im0,[]), title('Original Image');
subplot(122), imshow(ThisReconstruction{1},[]), title('Reconstructed Image');


% 3) Now let's synthesize a new image that should have the same parameters
% as the original one. Do not forget to include some parameters (max number
% of iterations, etc).
options.MaxIter = 20;
options.Size = size(im0);
[STex, Evolution] = TexSynthesizer(SPyr,TexFeats,options);


% Display the original and final result
figure, 
subplot(121), imshow(im0,[]); title('Original');
subplot(122), imshow(STex,[]); title('Synthetic');

% Display the evolution of the texture synthesis
figure, montage(Evolution./255)

% Or make a gif with the evolution
makeGIF(permute(Evolution./255,[1 2 4 3]),'myGIF.gif');

% Create a movie with the evolution of the texture synthesis
mov = immovie(repmat(Evolution./255,[1 1 3 1]));
implay(mov)


%% We can also create a synthesized texture AROUND the original one, by 
%  applying a mask. This process is equivalent to a texture extrapolation.
[M,N,~] = size(SPyr.Original);
Mask = zeros(M,N);
Mask(floor(M/4):M-floor(M/4),floor(N/4):N-floor(N/4)) = 1;
Mask = imgaussfilt(Mask,20);
figure, imshow(SPyr.Original.*Mask,[]); title('Image to be extrapolated');

[STex, Evolution] = TexSynthesizer(SPyr,TexFeats,options,Mask);

% Create a movie with the evolution of the texture synthesis
mov = immovie(repmat(Evolution./255,[1 1 3 1]));
implay(mov)