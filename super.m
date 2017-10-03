
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

function [ZInit, ZHR, ZInit_L, ZHR_L] = super( ...
        ZLR, ...
        factor, ...
        yResSubLR, xResSubLR, ...
        overlapLR, ...
        patRad, intSigma, dispMax, ...
        lambda0, lambda1, lambda2, ...
        innerCycles, outerCycles, ...
        warpMode, ...
        guessFlag, ...
        poolSize, ...
        alpha ...
        )
% SUPER super resolves the input light field. It divides the Low Resolution
% input light field into different sub light fields and super-resolves their
% luminances via GBSUPER. The super-resolved luminance of the different
% sub light fields are merged to produce the lumianance of the complete
% super-resolved light field. The chrominances of the input Low Resolution
% light field are simpy up-sampled via bilinear interpolation.
%
% INPUT:
% ZLR - the low resolution light field to super-resolve (double [0,1]).
% factor - the super-resolution factor (integer).
% yResSubLR - the vertical spatial resolution of a Low Resolution sub light field.
% xResSubLR - the horizontal spatial resolution of a Low Resolution sub light field.
% overlap - the pixel overlap between sub light fields.
% patRad - the size of the (patRad x patRad) patch used in the graph construction.
% intSigma - the standard deviation value used in the graph construction.
% dispMax - the maximum assumed disparity value.
% lambda0 - the weight of the data fidelity term.
% lambda1 - the weight of the warping term.
% lambda2 - the weight of the graph-based regularizer.
% innerCycles - the maximum number of iteration to be performed by the solver.
% outerCycles - the number of iterations of the GB super-resolution algorithm.
% warpMode - 0 for the Direct warping matrix construction, and 1 for the
%            SQuare-constraint matrix construction.
% guessFlag - 0 to initialize the solver with a zero light field, and 1 to
%             initialize the solver with ZInit.
% poolSize - the number of threads activated by PARFOR.
% alpha - the Tukey window parameter (the window is used in the merging phase).
%
% OUTPUT:
% ZInit - the RGB bilinear up-sampled version of ZLR (uint8 [0,255]).
% ZHR - the RGB Graph-Based super-resolved version of ZLR (uint8 [0,255]).
% ZInit_L - the luminance of ZInit (double [0,1]).
% ZHR_L - the luminance of ZHR (double [0,1]).
%
% Note that ZInit_L and ZHR_L are the light field luminances obtained from
% the up-sampling and super-resolution of ZLR, respectively. They are NOT
% extracted from ZInit and ZHR, as these two are in the uint8 data type,
% therefore they have undergone a quantization step.
%
% Note that this function assumes that the views in Z have dimensions which
% are multiple of the super-resolution factor.

% ==== Light field dimensions =============================================

% Angular resolution.
vRes = size(ZLR, 1);
hRes = size(ZLR, 2);
M = vRes * hRes;

% Spatial resolution (LR).
yResLR = size(ZLR{1, 1}, 1);
xResLR = size(ZLR{1, 1}, 2);

% Spatial resolution (HR).
yResHR = factor * yResLR;
xResHR = factor * xResLR;

% ==== RGB to LC1C2 =======================================================

% Turn ZLR from RGB to LC1C2.
[aux, lmin, lmax, c1min, c1max, c2min, c2max] = lc1c2(ZLR, 'forward');

% Separate the channels L, C1, and C2 ...

aux = reshape(lf2col(aux), [], 3);

ZLR_L = col2lf( ...
    aux(:, 1), ...
    vRes, hRes, ...
    yResLR, xResLR, 1 ...
    );

ZLR_C1 = col2lf( ...
    aux(:, 2), ...
    vRes, hRes, ...
    yResLR, xResLR, 1 ...
    );

ZLR_C2 = col2lf( ...
    aux(:, 3), ...
    vRes, hRes, ...
    yResLR, xResLR, 1 ...
    );

% ==== Compute an initial estimate of the HR light field ==================

% Perform interpolation of each view, separately.
aux = zeros(yResHR, xResHR);
ZInit_L = cell(vRes, hRes);
ZInit_C1 = cell(vRes, hRes);
ZInit_C2 = cell(vRes, hRes);
for u = 1:1:M
    
    aux(:, :) = upsample(ZLR_L{u}, factor, 'linear');
    aux(aux < lmin) = lmin;
    aux(aux > lmax) = lmax;
    ZInit_L{u} = aux;
    
    aux(:, :) = upsample(ZLR_C1{u}, factor, 'linear');
    aux(aux < c1min) = c1min;
    aux(aux > c1max) = c1max;
    ZInit_C1{u} = aux;
    
    aux(:, :) = upsample(ZLR_C2{u}, factor, 'linear');
    aux(aux < c2min) = c2min;
    aux(aux > c2max) = c2max;
    ZInit_C2{u} = aux;
    
end

% ==== Split ZLR_L into sub light fields ==================================

% Sub-light-field spatial resolution (HR).
yResSubHR = factor * yResSubLR;
xResSubHR = factor * xResSubLR;
overlapHR = factor * overlapLR;

% Split ZLR_L.
[zSubLR, ~, ~] = split(ZLR_L, yResSubLR, xResSubLR, overlapLR);

% Split the initial estimate ZInit_L.
[zSubInit, yCoord, xCoord] = ...
    split(ZInit_L, yResSubHR, xResSubHR, overlapHR);

% Number of sub-light-fields.
subNum = size(zSubLR, 2);

fprintf('Number of SUB-LIGHT-FIELDS: %d\n', subNum);

% ==== Sub-light-field blur and decimation matrices =======================

% Blur.
auxB = blurmat(yResSubHR, xResSubHR, factor);
B = kron(speye(M, M), auxB);

% Decimation.
auxD = decimat(yResSubHR, xResSubHR, factor);
D = kron(speye(M, M), auxD);

% ==== Super-resolve ZLR_L ================================================

% Super-resolve each LR sub-light-field, separately.
aux = cell(subNum, 1);
if poolSize > 1
    
    % PARALLEL COMPUTING ...
    
    % Initialize the pool.
    p = parpool(poolSize);
    
    parfor k = 1:1:subNum
        
        aux{k} = ...
            gbsuper( ...
            col2lf(zSubLR(:, k), vRes, hRes, yResSubLR, xResSubLR, 1), ...
            col2lf(zSubInit(:, k), vRes, hRes, yResSubHR, xResSubHR, 1), ...
            B, D, ...
            patRad, intSigma, dispMax, ...
            lambda0, lambda1, lambda2, ...
            lmin, lmax, ...
            innerCycles, outerCycles, ...
            warpMode, ...
            guessFlag);
    
    end
    
    % Terminate the pool.
    delete(p);
    
else
    
    % SEQUENTIAL COMPUTING ...
    
    for k = 1:1:subNum
        
        aux{k} = ...
            gbsuper( ...
            col2lf(zSubLR(:, k), vRes, hRes, yResSubLR, xResSubLR, 1), ...
            col2lf(zSubInit(:, k), vRes, hRes, yResSubHR, xResSubHR, 1), ...
            B, D, ...
            patRad, intSigma, dispMax, ...
            lambda0, lambda1, lambda2, ...
            lmin, lmax, ...
            innerCycles, outerCycles, ...
            warpMode, ...
            guessFlag);
    
    end
    
end

% Vectorize all the super-resolved sub-light-fields.
zSubHR = zeros(size(zSubInit, 1), subNum);
for k = 1:1:subNum
    
    zSubHR(:, k) = lf2col(aux{k});
    
end

% Merge the super-resolved sub-light-fields.
ZHR_L = merge( ...
    zSubHR, ...
    vRes, hRes, ...
    yResHR, xResHR, ...
    yResSubHR, xResSubHR, ...
    yCoord, xCoord, ...
    alpha ...
    );

% ==== Upsamples ZLR_C1 and ZLR_C2 ====================================

% Perform interpolation of each view, separately.
aux = zeros(yResHR, xResHR);
ZHR_C1 = cell(vRes, hRes);
ZHR_C2 = cell(vRes, hRes);
for u = 1:1:M
    
    aux(:, :) = upsample(ZLR_C1{u}, factor, 'cubic');
    aux(aux < c1min) = c1min;
    aux(aux > c1max) = c1max;
    ZHR_C1{u} = aux;
    
    aux(:, :) = upsample(ZLR_C2{u}, factor, 'cubic');
    aux(aux < c2min) = c2min;
    aux(aux > c2max) = c2max;
    ZHR_C2{u} = aux;
    
end

% ==== LC1C2 to RGB =======================================================

% Process ZInit ...

% Organize the channels L, C1, and C2 of ZInit into a single light field.
ZInit = col2lf( ...
    cat(1, lf2col(ZInit_L), lf2col(ZInit_C1), lf2col(ZInit_C2)), ...
    vRes, hRes, ...
    yResHR, xResHR, 3 ...
    );

% Converts ZInit from LC1C2 to RGB.
ZInit = lc1c2(ZInit, 'backward');

% Process ZHR ...

% Organize channel L, C1, and C2 of ZHR into a single light field.
ZHR = col2lf( ...
    cat(1, lf2col(ZHR_L), lf2col(ZHR_C1), lf2col(ZHR_C2)), ...
    vRes, hRes, ...
    yResHR, xResHR, 3 ...
    );

% Convert ZHR from LC1C2 to RGB.
ZHR = lc1c2(ZHR, 'backward');

end

