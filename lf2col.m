
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

function z = lf2col(Z)
% LF2COL vectorizes the input light field Z. The 2D cell array Z is scanned
% in column major order, and the same for its views Z{t,s}.
% If the light field views have more than one channel, then all the light
% fields (one for each channel) are vectorized separately, and then the
% vectorized light fields are stacked together.
%
% INPUT:
% Z - a light field.
%
% OUTPUT:
% z - the vectorized light field.

% =========================================================================

aux = cell2mat(Z(:)');
z = aux(:);

end

