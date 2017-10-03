
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

function P = psnrfr(A, B, peak, frame)
% PSNRFR computes the PSNR of image A with reference to image B, but it
% removes a frame around the image before the computation.
%
% INPUT:
% A - the image to evaluate,
% B - the reference image,
% peak - the maximum pixel value (e.g., 0 and 255),
% frame - the number of pixels to remove from the border.
%
% OUTPUT:
% P - the PSNR.

% =========================================================================

aux1 = A((1 + frame):(end - frame), (1 + frame):(end - frame), :);
aux2 = B((1 + frame):(end - frame), (1 + frame):(end - frame), :);

P = psnr(double(aux1), double(aux2), peak);

end

