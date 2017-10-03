
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

function P = psnrlf(Z1, Z2, peak, frame)
% PSNRLF computes the PSNR of the light field V1 with reference to the light field V2.
%
% INPUT:
% Z1 - the light field to evaluate,
% Z2 - the reference light field,
% peak - the maximum pixel value (e.g., 0 and 255),
% frame - the number of pixels to remove from the border.
%
% OUTPUT:
% P - a 2D matrix with P(i,j) the PSNR of the view V1{i,j}.

% =========================================================================

% Angular resolution.
vRes = size(Z1, 1);
hRes = size(Z1, 2);

P = zeros(vRes, hRes);
for s = 1:1:hRes
    for t = 1:1:vRes
        
        % Remove the border.
        aux1 = Z1{t, s}((1 + frame):(end - frame), (1 + frame):(end - frame), :);
        aux2 = Z2{t, s}((1 + frame):(end - frame), (1 + frame):(end - frame), :);
        
        % Compute the PSNR.
        P(t, s) = psnr(double(aux1), double(aux2), peak);
    
    end
end

end

