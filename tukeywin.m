
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

function W = tukeywin(N, alpha)
% TUKEYWIN builds a 1D Tukey window.
%
% INPUT:
% N - the window length.
% alpha - the window parameter.
%
% OUTPUT:
% W - the window.

% =========================================================================

n = 0:1:(N-1);
W = zeros(N, 1);

mask = (n < ceil(alpha * (N-1) / 2));
W(mask) = (1/2) * ( 1 + cos( pi * ( ( (2 * n(mask)) / (alpha * (N-1)) ) - 1 ) ) );

mask = (n >= ceil(alpha * (N-1) / 2)) & (n <= floor((N-1) * (1 - (alpha/2))));
W(mask) = 1;

mask = (n > floor((N-1) * (1 - (alpha/2))));
W(mask) = (1/2) * ( 1 + cos( pi * ( ( (2 * n(mask)) / (alpha * (N-1)) ) - (2/alpha) + 1 ) ) );

end

