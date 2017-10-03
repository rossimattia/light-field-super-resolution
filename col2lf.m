
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

function Z = col2lf(z, vRes, hRes, yRes, xRes, channels)
% COL2LF reverts the operation performed by function LF2COL.
%
% INPUT:
% z - the vectorized light field.
% vRes - the light field vertical angular resolution.
% hRes - the light field horizontal angular resolution.
% yRes - the light field vertical spatial resolution.
% xRes - the light field horizontal spatial resolution.
% channels - the number of channels in each view.
%
% OUTPUT:
% Z - the input light field stored as a (vRes x hRes) cell array of
%     (yRes x xRes x channels) views.

% =========================================================================

% Number of views.
M = vRes * hRes;

% Perform the reshaping.
aux = reshape(z, yRes, (length(z) / channels) / yRes, channels);
if (channels > 1)
    
    aux = mat2cell(aux, yRes, ones(M, 1) * xRes, channels);
    
else
    
    aux = mat2cell(aux, yRes, ones(M, 1) * xRes);
    
end
Z = reshape(aux, vRes, hRes);

end

