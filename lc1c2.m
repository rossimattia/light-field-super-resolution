
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

function [ZB, lmin, lmax, c1min, c1max, c2min, c2max] = lc1c2(ZA, direction)
% LC1C2 converts the input light field from RGB to LC1C2 if direction is
% 'forward', and viceversa if direction is 'backward'.
%
% INPUT:
% ZA - the input light field.
% direction - a string specifing the conversion type.
%
% OUTPUT:
% ZB - the transformed light field.
% lmin - the minimum value of the L channel for R, G, and B, in [0,1].
% lmax - the maximun value of the L channel for R, G, and B, in [0,1].
% c1min - the minimum value of the C1 channel for R, G, and B, in [0,1].
% c1max - the maximun value of the C1 channel for R, G, and B, in [0,1].
% c2min - the minimum value of the C2 channel for R, G, and B, in [0,1].
% c2max - the maximun value of the C2 channel for R, G, and B, in [0,1].

% ==== Trasformation matrices =============================================

% From RGB to LC1C2.
MF = ...
    [ ...
     1/4        1/2         1/4; ...
    -1/4        1/2        -1/4; ...
    -1/4        0.0         1/4 ...
    ];

% From LC1C2 to RGB.
MB = ...
    [ ...
    1.0        -1.0        -2.0; ...
    1.0         1.0         0.0; ...
    1.0        -1.0         2.0 ...
    ];

% Set the boundary values.
lmin = 0;
lmax = 1;
c1min = -1/2;
c1max = 1/2;
c2min = -1/4;
c2max = 1/4;

% ==== Light field dimensions =============================================    

% Angular resolution.
vRes = size(ZA, 1);
hRes = size(ZA, 2);

% Spatial resolution.
yRes = size(ZA{1, 1}, 1);
xRes = size(ZA{1, 1}, 2);

% Check that the input views have 3 channels.
if (size(ZA{1, 1}, 3) ~= 3)
    
    error('Input views must have 3 channels !!!\n\n');
    
end

% Select the correct transformation matrix.
switch direction
    
    case 'forward'
        
        M = MF;
        
    case 'backward'
        
        M = MB;
        
    otherwise
        
        error('Bad input arguments !!!\n\n');
        
end

% Perform the transformation ...

% Vectorize ZA.
zA = lf2col(ZA);

% Store the 3 vectorized light fields as the 3 columns of aux.
aux = reshape(zA, [], 3);

% Transform the 3 columns of aux.
aux = aux * (M');
% Note that the previous line is equivalent to:
% aux = (M * (aux'))';

% Arrange the tranformed (and vectorized) light field into a 2D cell array of views.
ZB = col2lf( ...
    aux(:), ...
    vRes, hRes, ...
    yRes, xRes, 3 ...
    );

end

