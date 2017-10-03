
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

function [D, newHeight, newWidth] = decimat(height, width, factor)
% DECIMAT builds a matrix that permits to implement image decimation as a
% matrix/vector multiplication.
%
% INPUT:
% height - the image height.
% width - the image width.
% factor - the decimation factor (sampling rate).
%
% OUTPUT:
% B - the decimation matrix.
%
% The decimation of an (height x width) image X can be carried out on the
% rows and the columns separately.
%
% Let Dc be the decimation matrix for the columns, and let Dr be the one
% for the rows. Then X can be decimated as follows:
% Y = (Dr*(Dc*X)')' = (Dr * X' * Dc')' = Dc * X * Dr'
%
% In vectorized form we have Y(:) = kron(Dr, Dc) * X(:), hence D = kron(Dr, Dc).

% ==== Compute the decimation matrix ======================================

% Compute Dc and Dr.
Dc = decimtx1D(height, factor);
Dr = decimtx1D(width, factor);

% Compute the decimation matrix D.
D = kron(Dr, Dc);
newHeight = size(Dc, 1);
newWidth= size(Dr, 1);

end

function mtx = decimtx1D(n, factor)

height = ceil(double(n) / factor);
width = n;

rows = zeros(height, 1);
cols = zeros(height, 1);

counter = 1;
for k = 1:1:height
    
    rows(k) = k;
    cols(k) = counter;
    
    counter = counter + factor;
    
end

mtx = sparse(rows, cols, ones(height, 1), height, width);

end

