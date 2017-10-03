
% =========================================================================
% =========================================================================
%
% Author:
% Mattia Rossi (rossi.mattia@gmail.com)
% Signal Processing Laboratory 4 (LTS4)
% Ecole Polytechnique Federale de Lausanne (Switzerland)
%
% =========================================================================
% =========================================================================

function [zSub, yCoord, xCoord] = split(Z, yResSub, xResSub, overlap)
% SPLIT divides the input light field into sub light fields.
%
% INPUT:
% Z - a light field.
% yResSub - the vertical spatial resolution of a sub light field.
% xResSub - the horizontal spatial resolution of a sub light field.
% overlap - the pixel overlap between sub light fields.
%
% OUTPUT:
% zSub - the sub light fields in vectorized form.
% yCoord - the vertical spatial coordinate of the top left pixel of each sub light field.
% xCoord - the horizontal spatial coordinate of the top left pixel of each sub light field.

% =========================================================================

% Light field Z angular resolution.
vRes = size(Z, 1);
hRes = size(Z, 2);
M = vRes * hRes;

% Light field Z spatial resolution.
yRes = size(Z{1, 1}, 1);
xRes = size(Z{1, 1}, 2);

% Compute the VERTICAL extraction coordinates for each view.
yStep = yResSub - overlap; 
y = [(1:yStep:(yRes - yResSub)), (yRes - yResSub + 1)];

% Compute the HORIZONTAL extraction coordinates for each view.
xStep = xResSub - overlap;
x = [(1:xStep:(xRes - xResSub)), (xRes - xResSub + 1)];

% Compute ALL the sub light field coordinates (in a view).
[X, Y] = meshgrid(x, y);
yCoord = Y(:);
xCoord = X(:);

% Compute the number of sub light fields.
subNum = length(yCoord);

% Allocate the space for the (vectorized) sub light fields.
zSub = zeros(yResSub * xResSub * M, subNum);

% Stack all the views, proceding in column major order.
stack = cell2mat( reshape(Z, [1, 1, M]) );

% Extract the sub light fields.
ptr = 1;
aux = zeros(yResSub, xResSub, M);
for k = 1:1:subNum
    
    aux(:, :, :) = stack( ...
        yCoord(k):(yCoord(k) + yResSub - 1), ...
        xCoord(k):(xCoord(k) + xResSub - 1), ...
        :);
    
    zSub(:, ptr) = aux(:);
    
    ptr = ptr + 1;
    
end

end

