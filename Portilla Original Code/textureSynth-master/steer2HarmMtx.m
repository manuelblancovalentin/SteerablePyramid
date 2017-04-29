% MTX = steer2HarmMtx(HARMONICS, ANGLES, REL_PHASES)
%
% Compute a steering matrix (maps a directional basis set onto the
% angular Fourier harmonics).  HARMONICS is a vector specifying the
% angular harmonics contained in the steerable basis/filters.  ANGLES 
% (optional) is a vector specifying the angular position of each filter.  
% REL_PHASES (optional, default = 'even') specifies whether the harmonics 
% are cosine or sine phase aligned about those positions.
% The result matrix is suitable for passing to the function STEER.

% Eero Simoncelli, 7/96.

function mtx = steer2HarmMtx(harmonics, angles, evenorodd)

%%=================================================================
%%% Optional Parameters:

if (exist('evenorodd') ~= 1)
  evenorodd = 'even';
end

% Make HARMONICS a row vector
harmonics = harmonics(:)';

numh = 2*size(harmonics,2) - any(harmonics == 0);

if (exist('angles') ~= 1)
  angles = pi * [0:numh-1]'/numh;
else
  angles = angles(:);
end

%%=================================================================

if isstr(evenorodd)
  if strcmp(evenorodd,'even')
    evenorodd = 0;
  elseif strcmp(evenorodd,'odd')
    evenorodd = 1;
  else
    error('EVEN_OR_ODD should be the string  EVEN or ODD');
  end
end

%% Compute inverse matrix, which maps Fourier components onto 
%% steerable basis.
imtx = permute(angles*harmonics,[1 3 2]);
if evenorodd
    imtx = reshape(cat(2,sin(imtx),-cos(imtx)),numel(angles),[],1);
else
    imtx = reshape(cat(2,cos(imtx),sin(imtx)),numel(angles),[],1);
end

  
r = rank(imtx);
if (( r ~= numh ) & ( r ~= size(angles,1) ))
  fprintf(2,'WARNING: matrix is not full rank');
end  

mtx = pinv(imtx);

