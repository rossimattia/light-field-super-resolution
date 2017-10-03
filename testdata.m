
% This script is meant to be run LOCALLY.
%
% It down-samples the datasets produced by MATLAB script "importdata", thus
% producing the LR input data for both my SR algorithm and COCOLIB-6.
%
% Note that the 2 folders "hci" and "stanford", produced by IMPORTDATA must
% be stored inside folder pathIn defined below.

clear;
close all;
clc;

% ==== Test parameters ====================================================

% HCI = 1 creates HCI LR data, HCI = 0 does not. The same for STF.
HCI = 1;
STF = 0;

% Down sampling factor.
factor = 2;

% Gaussian noise standard deviation.
sigmaNoise = 0.0;

% Initializes the seed for random number generation.
rng(0);

% Input/output paths.
pathIn = '/Users/mattiarossi/Documents/MATLAB/LFdatasets/';
pathOut = '/Users/mattiarossi/Documents/MATLAB/superdata/';

fprintf('SUB SAMPLING FACTOR IS %d\n', factor);

% ==== HCI dataset ========================================================

% Light field names.
if HCI
    
    fprintf('\nHCI dataset:\n');
    
    list = { ...
        'buddha', ...
        'buddha2', ...
        'couple' ...
        'cube', ...
        'horses', ...
        'maria' ...
        'medieval', ...
        'mona', ...
        'papillon', ...
        'pyramide', ...
        'statue', ...
        'stillLife' ...
        };
    
else
    
    list = {};
    
end

% Full input path.
pathInFull = [pathIn, 'hci/'];

% Full output path.
pathOutFull = [pathOut, 'hci', num2str(factor), '/'];
mkdir(pathOutFull);

% Down-samples each light field.
for k = 1:1:length(list)
    
    % Current light field name.
    name = list{k};
    
    fprintf('Sub sampling the light field %s.\n', name);
    
    % Loads the light field.
    aux = load([pathInFull, name, '.mat']);
    Z = aux.Z;
    clear('aux');
    
    % Sub samples the light field and adds random noise.
    [ZLR, ZHR] = high2low(Z, factor, sigmaNoise);
    
    % Saves the HR and LR light fields.
    save([pathOutFull, name, '.mat'], 'ZHR', 'ZLR');
    
    % Frees memory.
    clear('ZLR', 'ZHR');
    
end

% ==== STANFORD dataset ===================================================

% Light field names.
if STF
    
    fprintf('\nSTANFORD dataset:\n');
    
    list = { ...
        'amethyst', ...
        'beans', ...
        'bracelet', ...
        'bulldozer', ...
        'bunny', ...
        'cardsS', ...
        'chess', ...
        'eucalyptus', ...
        'knights', ...
        'treasure', ...
        'truck' ...
        };
    
else
    
    list = {};
    
end

% Full input path.
pathInFull = [pathIn, 'stanford/'];

% Full output path.
pathOutFull = [pathOut, 'stanford', num2str(factor), '/'];
mkdir(pathOutFull);

% Down-samples each light field.
for k = 1:1:length(list)
    
    % Current light field name.
    name = list{k};
    
    fprintf('Sub sampling the light field %s.\n', name);
    
    % Loads the light field.
    aux = load([pathInFull, name, '.mat']);
    Z = aux.Z;
    clear('aux');
    
    % Sub samples the light field.
    [ZLR, ZHR] = high2low(Z, factor, sigmaNoise);
    
    % Saves the HR and LR light fields.
    save([pathOutFull, name, '.mat'], 'ZHR', 'ZLR');
    
    % Frees memory.
    clear('ZLR', 'ZHR');
    
end

