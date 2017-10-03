
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

function L = stdlap(W)
% STDLAP computes the combinatorial Laplacian associated to the input similarity matrix.
% The Laplacian is defined as L = D-W, with D the graph degree matrix.
% D is diagonal with D(i,i) equal to sum(W(i, :)).
%
% INPUT:
% W - a similarity matrix.
%
% OUTPUT:
% L - the Laplacian matrix.

% =========================================================================

% Compute the degree matrix "D".
diagonal = sum(W, 2);
n = length(diagonal);
D = spdiags(diagonal, 0, n, n);

% Compute the Laplacian.
L = D - W;

end

