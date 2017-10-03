
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

function B = upsample(A, factor, method)
% UPSAMPLE up-samples the input image A by the value factor.
% The interpolation assumes that each input pixel has been obtained as the
% average of a (factor x factor) patch of pixels.
%
% INPUT:
% A - a gray scale image (double [0,1]).
% factor - the upsampling factor (integer),
% method - a string specifing the interpolation method to be used.
%
% OUTPUT:
% B - the up-sampled image.

% ==== Check input parameters =============================================

if (round(factor) - factor) ~= 0
    
    error('The upsampling factor must be integer !!!\n\n');
    
end

% ==== Dimensions =========================================================

% Input image resolution.
yResLR = size(A, 1);
xResLR = size(A, 2);

% Output image resolution.
yResHR = yResLR * factor;
xResHR = xResLR * factor;

% ==== Perform the interpolation ==========================================

% Compute the query coordinates.
[Xq, Yq] = meshgrid(1:1:xResHR, 1:1:yResHR);

% Compute the coordinates of the available samples.
firstX = (1 + factor) / 2;
firstY = firstX;
lastX = firstX + (factor * (xResLR - 1));
lastY = firstY + (factor * (yResLR - 1));
[X, Y] = meshgrid( ...
    (firstX - factor):factor:(lastX + factor), ...
    (firstY - factor):factor:(lastY + factor) ...
    );

% Create a 1 pixel frame around A.
aux = padarray(A, [1, 1], 'symmetric', 'both');

% Perform the interpolation.
B = interp2(X, Y, aux, Xq, Yq, method, 0);

end

