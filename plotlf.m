
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

function [] = plotlf(varargin)
% PLOTLF arranges the input light field view in a single image and plots it.
% In addition, it can plot (over the light field) the weights of the edges
% from a given pixel (t,s,y,x) to its neighbors.
%
% INPUT:
% varargin{1} - a light field, either gray scale or RGB.
% varargin{2} - the maximum intensity value in the light field.
% varargin{3} - an adjacency matrix of a graph on the light field.
% varargin{4} - the vertical angular coordinate of the pixel of interest.
% varargin{5} - the horizontal angular coordinate of the pixel of interest.
% varargin{6} - the vertical spatial coordinate of the pixel of interest.
% varargin{7} - the horizontal spatial coordinate of the pixel of interest.
%
% The input light field can be in any format (e.g., double, uint8, uint16).
% It will be converted in double [0,1] and plotted.

% ==== Constants ==========================================================

LEVELS = 2000;
SPACE = 1;

% ==== Check the number of input arguments ================================

if (nargin ~= 2) && (nargin ~= 7)
    
    error('Wrong input argument number !!!\n\n');
    
end

% ==== Read the input light field =========================================

% Read the light field and the peak value.
Z = varargin{1};
peak = varargin{2};

% Angular resolution.
vRes = size(Z, 1);
hRes = size(Z, 2);

% Spatial resolution.
yRes = size(Z{1, 1}, 1);
xRes = size(Z{1, 1}, 2);
N = yRes * xRes;

% Channels number (grays scale or RGB).
channels = size(Z{1, 1}, 3);

% ==== Processe the light field ===========================================

% Scale Z to [0,1] and extend it to RGB (if Z is gray scale).
Zrgb = cell(vRes, hRes);
if (channels == 1) || (channels == 3)
    
    Zrgb(:, :, :) = col2lf( ...
        repmat(double(lf2col(Z)) / peak, [4 - channels, 1]), ...
        vRes, hRes, ...
        yRes, xRes, 3 ...
        );
else
    
    error('Input light field in neither gray scale nor RGB !!!\n\m');

end

% ==== Processe the adjacency matrix ======================================

% If an adjacency matrix is provided too ...
if (nargin > 2)
    
    % Read the adjacency matrix.
    W = varargin{3};
    
    % Read the nodes whose connections we want to plot.
    t = varargin{4};
    s = varargin{5};
    y = varargin{6};
    x = varargin{7};
    
    % Compute the node linear index in the vectorised light field.
    k = ...
        N * ( sub2ind([vRes, hRes], t, s) - 1 ) + ...
        sub2ind([yRes, xRes], y, x);
    
    % Extract the "k"-th row of "W".
    wei = full(W(k, :))';
    
    % Define a colormap.
    map = zeros(LEVELS, 3);
    map(:, 1) = linspace(0, 1, LEVELS);
    
    % Turn wei values into colormap indexes.
    auxMax = max(wei);
    if auxMax > 0
        
        wei = round( ...
            ((LEVELS - 1) * (wei ./ auxMax)) + 1 ...
            );
        
    else
        
        % Avoid division by zero in the case auxMax is zero.
        wei = ones(length(wei));
        
    end
    
    % Convert wei from indexes of map rows to RGB values.
    wei = ind2rgb(wei, map);
    
    % Arrange wei RGB values into an RGB light field.
    Wrgb = col2lf(wei(:), vRes, hRes, yRes, xRes, 3);
    
    % Merges Zrgb and Wrgb ...
    
    % Vectorizes Zrgb and Wrgb.
    auxZ = lf2col(Zrgb);
    auxW = lf2col(Wrgb);
    
    % Detect those pixels that are connected (i.e., NON zero weight) to pixel k.
    aux = reshape(auxW, [], 3);
    aux = sum(aux, 2);
    mask = (aux > 0);
    
    % In Zrgb, replace the connected pixels with the corresponding weights of "Wrgb".
    mask = repmat(mask, [3, 1]);    % mask extension to the 3 color channels.
    auxZ(mask) = auxW(mask);
    
    % Rearrange auxZ as a light field.
    Zrgb(:, :) = col2lf(auxZ, vRes, hRes, yRes, xRes, 3);
    
end

% ==== Stitche together all the views in Zrgb  ============================

% Allocate the background.
I = zeros(vRes * (yRes + SPACE) + SPACE, hRes * (xRes + SPACE) + SPACE, 3);

% Place the views on the background I.
vOff = 1;
hOff = 1;
for s = 1:1:hRes
    for t = 1:1:vRes
        
        I((vOff + 1):(vOff + yRes), (hOff + 1):(hOff + xRes), :) = Zrgb{t, s};
        
        vOff = vOff + yRes + SPACE;
        
    end
    
    vOff = SPACE;
    hOff = hOff + xRes + SPACE;
    
end

% ==== Plot the stitched viwes ============================================

imshow(I);

end

