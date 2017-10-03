
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

function B = blurmat(height, width, factor)
% BLURMAT builds a matrix that permits to implement image blurring as a
% matrix/vector multiplication. The blurring consists in a box kernel.
%
% INPUT:
% height - the image height.
% width - the image width.
% factor - the size of the (factor x factor) blurring kernel.
%
% OUTPUT:
% B - the blur matrix.
%
% The blurring of an image is defined as the convolution of the image with
% a 2D box kernel. Since 2D box kernels are separable, convolution of a
% (height x width) image X with a 2D box kernel can be implemented by first
% filtering the columns of X with a 1D box kernel, and then its rows with
% the same kernel.
%
% Let Bc be the filtering matrix for the columns, and let Br be the one for
% the rows. Then X can be filtered as follows:
% Y = (Br*(Bc*X)')' = (Br * X' * Bc')' = Bc * X * Br'
%
% In vectorized form we have Y(:) = kron(Br, Bc) * X(:), hence B = kron(Br, Bc).
%
% Note that the matrices Bc and Br may be different, despite implementing
% the same kernel, as image X may not be square.

% ==== Compute the blur matrix ============================================

% Compute Bc and Br.
Bc = avgmtx(height, factor);
Br = avgmtx(width, factor);

% Compute the blur matrix B.
B = kron(Br, Bc);

end

function mtx = avgmtx(n, factor)

kernel = zeros((2 * (factor - 1)) + 1, 1);
kernel(factor:end) = 1 / factor;

diags = -(factor - 1):1:(factor - 1);

D = repmat(kernel', [n, 1]);
mtx = spdiags(D, diags, n, n);

end

