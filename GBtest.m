
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
%
% This script permits to test the Graph-Based Light Field Super-Resolution
% algorithm, presented in the papers below, on the HCI and Stanford light
% field datasets.
%
% Mattia Rossi and Pascal Frossard
% "Graph-Based Light Field Super-Resolution"
% in 19th IEEE International Workshop on Multimedia Signal Processing, 2017.
%
% Mattia Rossi and Pascal Frossard
%"Light Field Super-Resolution Via graph-Based Regularization"
% https://arxiv.org/pdf/1701.02141.pdf
%
% The "TEST PARAMETERS" section allows the user to play with the algorith
% parameters, such as the super-resolution factor, the approach to the
% warping matrix construction, etc.
%
% =========================================================================
% =========================================================================

clear;
close all;
clc;

%  ==== TEST PARAMETERS ===================================================

% "HCI = 1" allows tests on HCI data, HCI = 0 does not. The same for STF.
HCI = 1;
STF = 1;

% Warping mode (DR or SQ).
warpMode = 'SQ';
%warpMode = 'DR';

% Set the input/ouput paths.
pathIn = '/Users/mattiarossi/Documents/MATLAB/superdata/';
pathOut = '';

% Set the PARFOR pool size.
poolSize = 1;

% Super resolution factor.
factor = 2;

% Desired angular resolution.
vRes = 5;
hRes = vRes;

% HCI angular cropping parameters.
tminHCI = 3;
tmaxHCI = tminHCI + vRes - 1;
sminHCI = 3;
smaxHCI = sminHCI + hRes - 1;

% STANFORD angular cropping parameters.
tminSTF = 7;
tmaxSTF = tminSTF + vRes - 1;
sminSTF = 7;
smaxSTF = sminSTF + hRes - 1;

% Spatial cropping parameters (at LR). This block allows the user to
% super-resolve only a small portion of the Low Resolution views. This
% could be useful when playing with the algorithm parameters.
%offSet = 80;
%yStartLR = 70;
%xStartLR = 100;
%yRangeLR = [yStartLR, yStartLR + offSet - 1];
%xRangeLR = [xStartLR, xStartLR + offSet - 1];
%yAux = ((yStartLR - 1) * factor);
%xAux = ((xStartLR - 1) * factor);
%yRangeHR = [yAux + 1, yAux + (offSet * factor)];
%xRangeHR = [xAux + 1, xAux + (offSet * factor)];
yRangeLR = [1, -1]; % No cropping.
xRangeLR = [1, -1]; % No cropping.
yRangeHR = [1, -1]; % No cropping.
xRangeHR = [1, -1]; % No cropping.

% Inter (views) graph parameters:
patRad = 3;             % Square patch radius.
intSigma = 0.7229;      % Weights decay.
dispMax = 6;            % Maximum disparity.

% Energy function multipliers.
lambda0 = 1.0;
lambda1 = 0.2;
lambda2 = 0.0055;

% Number of inner and outer iterations. 
innerCycles = 200;
outerCycles = 2;

% Activation flag for the initial guess.
guessFlag = 1;

% Spatial resolution for each LR sub light field.
yResSubLR = 100;
xResSubLR = 100; 

% Pixel overlap between two adjacent LR sub light fields.
overlapLR = ceil(yResSubLR / 2);

% Merging window parameter.
alpha = 0.15;

list = { ...
    {'buddha',          'hci'}, ...
    {'buddha2',         'hci'}, ...
    {'couple',          'hci'}, ...
    {'cube',            'hci'}, ...
    {'horses',          'hci'}, ...
    {'maria',           'hci'}, ...
    {'medieval',        'hci'}, ...
    {'mona',            'hci'}, ...
    {'papillon',        'hci'}, ...
    {'pyramide',        'hci'}, ...
    {'statue',          'hci'}, ...
    {'stillLife',       'hci'}, ...
    {'amethyst',        'stanford'}, ...
    {'beans',           'stanford'}, ...
    {'bracelet',        'stanford'}, ...
    {'bulldozer',       'stanford'}, ...
    {'bunny',           'stanford'}, ...
    {'cardsS',          'stanford'}, ...
    {'chess',           'stanford'}, ...
    {'eucalyptus',      'stanford'}, ...
    {'knights',         'stanford'}, ...
    {'treasure',        'stanford'}, ...
    {'truck',           'stanford'} ...
    };

% ==== SUPER-RESOLUTION ===================================================

% Super resolve each light field.
for k = 1:1:length(list)
    
    % Light field name and dataset.
    name = list{k}{1};
    dataset = list{k}{2};
    
    switch dataset
        
        case 'hci'
            
            if (HCI > 0)
                
                tRange = [tminHCI, tmaxHCI];
                sRange = [sminHCI, smaxHCI];
                
            else
                
                continue;
                
            end
            
        case 'stanford'
            
            if (STF > 0)
                
                tRange = [tminSTF, tmaxSTF];
                sRange = [sminSTF, smaxSTF];
                
            else
                
                continue;
                
            end
            
        otherwise
            
            error('Dataset not supported !!!\n\n');
            
    end
    
    % Full input path.
    pathInFull = [pathIn, dataset, num2str(factor), '/'];
    
    fprintf('Processing the light field %s.\n', upper(name));
    
    % Load the test light field.
    aux = load([pathInFull, name, '.mat']);
    ZHR = aux.ZHR;
    ZLR = aux.ZLR;
    clear('aux');
    
    % Crop the angular resolution of ZHR and ZLR.
    ZHR = crop(ZHR, tRange, sRange, yRangeHR, xRangeHR);
    ZLR = crop(ZLR, tRange, sRange, yRangeLR, xRangeLR);
    
    % Spatial resolution HR.
    yResHR = size(ZHR{1, 1}, 1);
    xResHR = size(ZHR{1, 1}, 2);
    
    % Spatial resolution LR.
    yResLR = size(ZLR{1, 1}, 1);
    xResLR = size(ZLR{1, 1}, 2);
    
    if ((yResSubLR > yResLR) || (xResSubLR > xResLR))
        
        error('The spatial resolution of a sub light field cannot exceed the one of the full light field !!! \n\n');
        
    end
    
    % Super-resolve the light field ZLR. The computation is performed with pixels in double [0,1].
    timeStart = tic;
    [ZInit, ZGB, ZInit_L, ZGB_L] = super( ...
        col2lf(double(lf2col(ZLR)) / 255, vRes, hRes, yResLR, xResLR, 3), ...   % ZLR to double [0,1].
        factor, ...
        yResSubLR, xResSubLR, ...
        overlapLR, ...
        patRad, intSigma, dispMax, ...
        lambda0, lambda1, lambda2, ...
        innerCycles, outerCycles, ...
        warpMode, ...
        guessFlag, ...
        poolSize, ...
        alpha);
    timeElapsed = toc(timeStart);
    
    % Note that ZInit, ZGB, ZInit_L, ZGB_L are double [0,1].
    
    % Quantize ZInit pixels to uint8 [0,255].
    ZInit = col2lf( ...
        uint8(round(lf2col(ZInit) * 255)), ...
        vRes, hRes, ...
        yResHR, xResHR, 3 ...
        );
    
    % Quantize ZGB pixels to uint8 [0,255].
    ZGB = col2lf( ...
        uint8(round(lf2col(ZGB) * 255)), ...
        vRes, hRes, ...
        yResHR, xResHR, 3 ...
        );
    
    % ==== Compute the PSNRs ==============================================
    
    % Extract the luminance of the original High Resolution light field.
    aux = double(lf2col(ZHR)) / 255;
    aux = col2lf(aux, size(ZHR, 1), size(ZHR, 2), size(ZHR{1, 1}, 1), size(ZHR{1, 1}, 2), 3);
    aux = lc1c2(aux, 'forward');
    aux = lf2col(aux);
    aux = aux(1:(size(ZHR, 1) * size(ZHR, 2) * size(ZHR{1, 1}, 1) * size(ZHR{1, 1}, 2)));
    aux = col2lf(aux, size(ZHR, 1), size(ZHR, 2), size(ZHR{1, 1}, 1), size(ZHR{1, 1}, 2), 1);
    
    % PSNR of the bilinearly interpolated light field.
    PInit = psnrlf(ZInit_L, aux, 1, 0);
    PInitMean = mean(PInit(:));
    
    % PSNR of the super-resolved light field.
    PGB = psnrlf(ZGB_L, aux, 1, 0);
    PGBMean = mean(PGB(:));
    
    % ==== Save the results ===============================================
    
    % Full output path.
    pathOutFull = [pathOut, dataset, '_sr', num2str(factor), '_', warpMode, '/'];
    mkdir(pathOutFull);
    
    % Store all the test parameters in a struct.
    metadata.warpMode = warpMode;
    metadata.pathInFull = pathInFull;
    metadata.pathOutFull = pathOutFull;
    metadata.poolSize = poolSize;
    metadata.factor = factor;
    metadata.name = name;
    metadata.dataset = dataset;
    metadata.tRange = tRange;
    metadata.sRange = sRange;
    metadata.yRangeLR = yRangeLR;
    metadata.xRangeLR = xRangeLR;
    metadata.yRangeHR = yRangeHR;
    metadata.xRangeHR = xRangeHR;
    metadata.patRad = patRad;
    metadata.intSigma = intSigma;
    metadata.dispMax = dispMax;
    metadata.lambda0 = lambda0;
    metadata.lambda1 = lambda1;
    metadata.lambda2 = lambda2;
    metadata.innerCycles = innerCycles;
    metadata.outerCycles = outerCycles;
    metadata.guessFlag = guessFlag;
    metadata.yResSubLR = yResSubLR;
    metadata.xResSubLR = xResSubLR;
    metadata.overlapLR = overlapLR;
    metadata.alpha = alpha;
    metadata.timeElapsed = timeElapsed;
    
    % Store the super-resolved light field.
    save( ...
        strcat(pathOutFull, name, '.mat'), ...
        'ZHR', 'ZLR', 'ZInit', 'ZGB', 'ZInit_L', 'ZGB_L', 'metadata' ...
        );
    
    % Plot the results.
    figure;
    
    subplot(2, 2, 1);
    plotlf(ZLR, 255);
    title('LR');
    
    subplot(2, 2, 2);
    plotlf(ZHR, 255);
    title('HR - Original');
    
    subplot(2, 2, 3);
    plotlf(ZInit, 255);
    str = sprintf('HR - Bilinear interp. %2.2f dB', PInitMean);
    title(str);
    
    subplot(2, 2, 4);
    plotlf(ZGB, 255);
    str = sprintf('HR - Graph-Based %2.2f dB', PGBMean);
    title(str);
    
    % Remove the big variables.
    clear( ...
        'ZHR', ...
        'ZLR',...
        'ZInit', 'ZInit_L', ...
        'ZEst', 'ZEst_L', ...
        'metadata' ...
        );
    
    fprintf('\n');
    
end

