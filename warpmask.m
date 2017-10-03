
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

function H = warpmask(W, B, D, vRes, hRes)
% WARPMASK computes the visibility masks for the warping. Let us assume that
% the matrix W corresponds to the j-th view. Then W contains all the warping
% matrices Wij (warping from the j-th view to the i-th one) stacked along
% the second dimension, for all the views i in the light field. WARPMASK
% computes a mask for each Wij and stores it in H{i}. The masks handle
% exclusively the dis/occlusions caused by the view borders.
%
% INPUT:
% W - a matrix containing the warping matrices (stacked along the second dimension)
%     from a given light field view to all the others.
% B - the blur matrix for the whole light field.
% D - the decimation matrix for the whole light field.
% vRes - the light field vertical angular resolution.
% hRes - the light field horizontal angular resolution.
%
% OUTPUT:
% H - a 1D cell array containing a mask matrix in each cell.

% =========================================================================

% Number of views.
M = vRes * hRes;

% Compute the warping mask (at Low Resolution) for all the warpings Wij.
% The computed masks are in vectorized form.
S = sum(W, 2);
S = spones(S);
S = D * B * S;
S(S < 1) = 0;

% Number of pixels at each LR view.
NLR = size(D, 1) / M;

% From the warping mask, create the corresponding diagonal matrices.
H = cell(M, 1);
for u = 1:1:M
    
    H{u} = spdiags( ...
        S((((u - 1) * NLR) + 1):(u * NLR)), ...
        0, ...
        NLR, NLR ...
        );

end

end

