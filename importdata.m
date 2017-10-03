
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

% This script reads each light field in the HCI and Stanford light field
% datasets, arranges them according to the conventions stated in the README
% file, and finally stores each of them in a separate ".mat" file.
% The ".mat" files are stored in the folders "pathOut/hci" and "pathOut/stanford",
% with pathOut defined below.

clear;
close all;
clc;

% =========================================================================

% Input path.
pathIn = '/Users/mattiarossi/Documents/Datasets/LFdatasets/';

% Output path.
pathOut = '/Users/mattiarossi/Documents/MATLAB/LFdatasets/';

% ==== HCI dataset ========================================================

fprintf('HCI dataset:\n');

% Light field list.
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

% Full input path.
pathInFull = [pathIn, 'hci/'];

% Full output path.
pathOutFull = [pathOut, 'hci/'];
mkdir(pathOutFull);

% Extract the light fields.
for k = 1:1:length(list)
    
    % Current light field name.
    name = list{k};
    fprintf('Processing the light field %s.\n', name);
    
    % Read the light field.
    lf = readhci([pathInFull, name, '/', 'lf.h5']);
    
    % Extract the views.
    Z = lf.view;
    
    % Save the light field views.
    save([pathOutFull, name, '.mat'], 'Z');
    
end

% ==== STANFORD dataset ===================================================

fprintf('\nSTANFORD dataset:\n');

% Light field list.
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

% Full input path.
pathInFull = [pathIn, 'stanford/'];

% Full output path.
pathOutFull = [pathOut, 'stanford/'];
mkdir(pathOutFull);

% Extract the light fields.
for k = 1:1:length(list)
    
    % Current light field name.
    name = list{k};
    fprintf('Processing the light field %s.\n', name);
    
    % Read the light field and extracts its views.
    Z = readstf([pathInFull, name, '/']);
    
    % Three light fields are stored differently: they need a permutation.
    if (strcmp(name, 'cardsL') || ...
            strcmp(name, 'cardsS') || ...
            strcmp(name, 'knights'))
        
        Z = fliplr(Z);
        
    end
    
    % Save the light field views.
    save([pathOutFull, name, '.mat'], 'Z');
    
end

