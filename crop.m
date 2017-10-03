
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

function ZCR = crop(Z, vRange, hRange, yRange, xRange)
% CROP reduces the angular and/or spatial resolution of the input light
% field by cropping the number of views (angular cropping) and/or the views
% area (spatial cropping).
%
% INPUT:
% Z - the light field to be cropped.
% vRange - a two element vector with the starting and ending vertical angular cropping indexes.
% hRange - a two element vector with the starting and ending horizontal angular cropping indexes.
% yRange - a two element vector with the starting and ending vertical spatial cropping indexes.
% xRange - a two element vector with the starting and ending horizontal spatial cropping indexes.
%
% OUTPUT:
% ZCR - the cropped light field.

% ==== Auxiliary variables ================================================

% Compute the (selection) vectors for angular cropping.
vSelec = vRange(1):1:vRange(2);
hSelec = hRange(1):1:hRange(2);

% Compute the (selection) vectors for spatial cropping.
ySelec = yRange(1):1:yRange(2);
xSelec = xRange(1):1:xRange(2);

% ==== Copy the input light field =========================================

ZCR = Z;

% ==== Perform the angular cropping =======================================

if ~(isempty(vSelec) || isempty(hSelec))
    
    ZCR = ZCR(vSelec, hSelec);
    
end

% ==== Perform the spatial cropping =======================================

if ~(isempty(ySelec) || isempty(xSelec))
    
    M = size(ZCR, 1) * size(ZCR, 2);
    
    for u = 1:1:M
        
        ZCR{u} = ZCR{u}(ySelec, xSelec, :);
    
    end
    
end

end

