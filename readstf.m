
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

function Z = readstf(path)
% READSTF reads the STANFORD light field specified in the input folder and
% it arranges it in a 2D cell array adopting the reference system described
% in the README.
%
% INPUT:
% path - the name of a folder containing the light field views.
%
% OUTPUT:
% Z - a light field in the format described in the README.

% =========================================================================

% Each STANFORD light field is a 16x16 grid of views.
vRes = 17;
hRes = 17;

% Parametrized name of a view.
str = [path, 'view_%02d_%02d.png'];

% Read the views paying attention to STANFORD reference system.
Z = cell(vRes, hRes);
for t = 1:1:vRes
    for s = 1:1:hRes
        
        name = sprintf(str, vRes - t, s - 1);
        Z{t, s} = imread(name);
    
    end
end

end

