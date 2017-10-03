
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

function [] = lf2grid(Z, path)
% LF2GRID receives a light field and saves it a set of PNG images named
% according to the reference system adopted in the HCI light field dataset.
%
% INPUT:
% Z - a light field.
% path - the destination path.

% =========================================================================

% Angular resolution.
vRes = size(Z, 1);
hRes = size(Z, 2);

% Create the output folder (in the case that it does not exist yet).
mkdir(path);

% Write each view to a PNG file, according to the HCI dataset convention.
for s = 1:1:hRes
    for t = 1:1:vRes
        
        name = sprintf([path, 'out_%02d_%02d.png'], s - 1, t - 1);
        imwrite(Z{vRes - t + 1, hRes - s + 1}, name);
    
    end
end

end

