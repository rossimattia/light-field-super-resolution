
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

function [ZLR, ZHR] = high2low(Z, factor, sigmaNoise)
% HIGH2LOW down-samples the light field views in Z by factor and adds
% Gaussian random noise with standard deviation sigmaNoise.
%
% INPUT:
% Z - the light field to down-sample.
% factor - the down-sampling factor.
% sigmaNoise - the Gaussian noise standard deviation.
%
% OUTPUT:
% ZLR - the low resolution light field.
% ZHR - the high resolution light field cropped such that its spatial
%       dimensions are multiple of the down-sampling factor.
%
% The output low resolution light field is to be considered the low resolution
% version of the OUTPUT high resolution light field, and NOT of the INPUT
% one, as only the output high resolution light field is properly cropped.

% ==== Check the input type ===============================================

if ~isa(lf2col(Z), 'uint8')
    
    error('Input light field must have uint8 views !!!\n\n');
    
end

% ==== Ligth field parameters =============================================

% Angular resolution.
vRes = size(Z, 1);
hRes = size(Z, 2);

% Spatial resolution.
yRes = size(Z{1, 1}, 1);
xRes = size(Z{1, 1}, 2);

% Channels number (gray scale or RGB).
channels = size(Z{1, 1}, 3);

% ==== Crop Z =============================================================

% Each view of Z is cropped at the right and bottom sides, such that their
% dimensions are multiples of factor.

% Compute the LR dimensions of the views.
yResLR = floor(yRes / factor);
xResLR = floor(xRes / factor);

% New spatial resolution.
yResHR = factor * yResLR;
xResHR = factor * xResLR;

% Perform the cropping.
ZHR = cell(vRes, hRes);
for s = 1:1:hRes
    for t = 1:1:vRes
        
        ZHR{t, s} = Z{t, s}(1:yResHR, 1:xResHR, :);
        
    end
end

% ==== Blur and Decimate "ZHR" ==========================================

% Blur matrix for a single view and channel.
B = blurmat(yResHR, xResHR, factor);

% Decimation matrix for a single view and channel.
D = decimat(yResHR, xResHR, factor);

% Blurring and decimation matrix for a single view and channel.
DB = D * B;

% Blurring and decimation matrix for a single view and ALL its channels.
DBch = kron(speye(channels), DB);

% ZLR will contain the Low Resolution (LR) views.
ZLR = cell(vRes, hRes);

% Number of pixels in each LR view.
n = yResLR * xResLR * channels;

% Blur, decimate, and add noise to each view, separately.
for s = 1:1:hRes
    for t = 1:1:vRes
        
        auxLR = (DBch * double(ZHR{t, s}(:))) + (sigmaNoise * 255 * randn(n, 1));
        auxLR = uint8(round(auxLR));
        
        ZLR{t, s} = reshape( ...
            auxLR, ...
            yResLR, ...
            xResLR, ...
            channels ...
            );
        
    end
end

end

