
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

function lf = readhci(filename)
% READHCI stores the data inside the HDF5 file at the input into a struct
% with the same fields. Views are arranged in a 2D cell array adopting the
% reference system described in the README.
%
% INPUT:
% filename - an HDF5 file name.
%
% OUTPUT:
% lf - a struct with:
% - the same fields of the HDF5 file,
% - a 2D cell array storing the light field,
% - a 2D cell array storing the light field depth maps.

% ==== Extract the light field parameters =================================

lf.yRes = h5readatt(filename, '/', 'yRes');
lf.xRes = h5readatt(filename, '/', 'xRes');
lf.vRes = h5readatt(filename, '/', 'vRes');
lf.hRes = h5readatt(filename, '/', 'hRes');
lf.channels = h5readatt(filename, '/', 'channels');
lf.vSampling = h5readatt(filename, '/', 'vSampling');
lf.hSampling = h5readatt(filename, '/', 'hSampling');
lf.focalLength = h5readatt(filename, '/', 'focalLength');
lf.dV = h5readatt(filename, '/', 'dV');
lf.dH = h5readatt(filename, '/', 'dH');
lf.shift = h5readatt(filename, '/', 'shift');

% ==== Check the input light field ========================================

if (lf.channels ~= 3)
    
    error('Input views must be RGB !!!');

end

% ==== Read the light field views =========================================

% Dataset '/LF' is a 5D array of size (channels x xRes x yRes x hRes x vRes).

% Read the light field views.
dataZ = h5read(filename, '/LF');

% Order dataLF as (yRes x xRes x channels x vRes x hRes).
dataZ = permute(dataZ, [3 2 1 5 4]);

% ==== Read the light field depth maps ====================================

% '/GT_DEPTH' is a 4D array of size (xRes x yRes x hRes x vRes).

% Read the light field depth maps.
dataW = h5read(filename, '/GT_DEPTH');

% Order dataW as (yRes x xRes x vRes x hRes).
dataW = permute(dataW, [2 1 4 3]);

% ==== Organize the views and depth maps into 2D cell arrays ==============

% Allocate a 2D cell array for the light field views.
lf.view = cell(lf.vRes, lf.hRes);

% Allocate a 2D cell array for the light field depth maps.
lf.depth = cell(lf.vRes, lf.hRes);

for s = 1:1:lf.hRes
    for t = 1:1:lf.vRes
        
        % Extract view (t,s) and its depth map.
        auxZ = dataZ(:, :, :, t, lf.hRes - s + 1);
        auxW = dataW(:, :, t, lf.hRes - s + 1);
        
        % Store the view and its depth map.
        lf.view{t, s} = auxZ;
        lf.depth{t, s} = auxW;
    
    end
    
end

end

