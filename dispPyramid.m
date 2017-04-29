function dispPyramid(ThisImage, ThisPyramid, TheseMarginals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function displays all the levels of the steerable pyramid in one 
% single image
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

if nargin < 3
    TheseMarginals = [];
end

ThisImage = mat2gray(ThisImage);

if size(ThisImage,3) ~= 3
    ThisImage = repmat(ThisImage,[1 1 3]);
end

if isempty(ThisPyramid), return; end

NLevels = numel(ThisPyramid);
NBands = size(ThisPyramid{1},3);

if isempty(TheseMarginals)
   for i=1:NLevels
        TheseMarginals{i} = imresize(ThisImage,1/2^(i-1)); 
   end
else
    for i=1:NLevels
        if size(TheseMarginals{i},3) ~= 3
            TheseMarginals{i} = repmat(TheseMarginals{i},[1 1 3]);
        end
        TheseMarginals{i} = mat2gray(real(TheseMarginals{i}));
    end
end

[M,N,K] = size(TheseMarginals{1});
% create the image first and then we will fill it in a loop
%ncols = floor(sqrt(2*round(NLevels/2)));
ncols = floor(sqrt(NBands));
nrows = NBands - ncols;

% add the row and col where the original image will be:
ncols = ncols + 1;
nrows = nrows + 1;

% Set the colors we are going to use as borders
cols = permute(label2rgb((1:NLevels)'),[1 3 2]);

% Initialize collage
Collage = zeros(ncols*M,nrows*N,K);

[M0,N0,~] = size(Collage);

for i=1:NLevels
   
    % original image in the right-bottom corner
    SubImage = TheseMarginals{i};
    [subM,subN,~] = size(SubImage);
    Collage(M0-subM+1:M0,N0-subN+1:N0,:) = SubImage;
    
    % from 1 to ncols-1, we plot on cols (left of subimage)
    % from ncols+1 : NBands we plot on rows (up of subimage)
    for j=1:ncols-1
        ist = M0-subM+1;
        ien = M0;
        jst = N0-(ncols-j+1)*subN + 1;
        jen = jst + subN - 1;
        Collage(ist:ien,jst:jen,:) = repmat(mat2gray(real(ThisPyramid{i}(:,:,j))),[1 1 K]);
    end
    for j=ncols:NBands
        ien = M0-(j-ncols+1)*subM + 1;
        ist = ien - subM + 1;
        jst = N0-subN+1;
        jen = N0;
        Collage(ist:ien,jst:jen,:) = repmat(mat2gray(real(ThisPyramid{i}(:,:,j))),[1 1 K]);
    end
    
    % add a rectangle so it's easier to visualize the collage
    mup = M0 - nrows*subM + 1;
    mdown = M0;
    mmiddle = mdown - subM;
    nleft = N0 - ncols*subN + 1;
    nright = N0;
    nmiddle = nright - subN;
    
    ThisPolygon = [nleft, mdown, nright, mdown, nright, mup,...
                   nmiddle, mup, nmiddle, mmiddle, nleft, mmiddle];
    
    Collage = insertShape(Collage,'polygon',ThisPolygon,...
                     'LineWidth',2,'Color',cols(i,:));
        

    M0 = M0-subM+1;
    N0 = N0-subN+1;
    
end


figure, 
imshow(Collage); axis normal;
title(sprintf('Steerable Pyramid Collage. #Levels = %i / #Bands = %i',NLevels,NBands));


end