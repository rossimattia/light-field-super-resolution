
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

function wei = nlm(A, B, patRad, intSigma, offset)
% NLM computes the similarity weights between the pixels of an image A and
% those of an image B. In graph terms, it computes the edge weights of a
% graph were edges go from the pixels of image A to those of image B. The
% weights are stored in a vector wei with three columns: the weight of the
% edge from pixel wei(k,1) of image A to the pixel wei(k,2) of image B is
% wei(k,3), with wei(k,1) and wei(k,2) linear coordinates.
%
% INPUT:
% A - a gray scale image (double [0,1]).
% B - a gray scale image (double [0,1]).
% patRad - the size of the (patRad x patRad) patch used in the similarity computation.
% intSigma - the standard deviation used in the similarity computation.
% offSet - a struct defining the searching window (see below).
%
% OUTPUT:
% wei - the similarity weights (aka the graph weights).
%
% Each pixel A(y,x) is compared with all the pixels B(y + offSet.Y(k),x + offSet.X(k))
% where k = 1, 2, ... , length(offSet.X). The comparison is carried out
% between the (patRad x patRad) patches centered at the two considered pixels.
%
% For a matter of computational efficiency, some rows of wei are computed
% but NOT valid. These rows are not valid because either wei(k,1) or wei(k,2)
% does not correspond to a real pixel coordinate in A or B, respectively.
% These rows are the only ones having wei(k,3) equal to ZERO, therefore
% they can be easily recognized.

% =========================================================================

tStart = tic;

% ==== Global constants ===================================================

QUIET = 1;

% ==== Extract the images dimensions ======================================

height = size(A, 1);
width = size(A, 2);
N = height * width;

% ==== Compute the weights ================================================

% Extract the comparison offsets.
Y = offset.Y;
X = offset.X;
offsetsNum = length(X);

% Compute the radius of the smallest window that contains the offsets.
% The window will be ((2 * winRad + 1) x (2 * winRad + 1)) pixels.
winRad = max( max(abs(Y)), max(abs(X)) );

% Pad images A and B.
Apad = padarray(A, [patRad, patRad], 0.5, 'both');
Bpad = padarray(B, [winRad + patRad, winRad + patRad], 0.5, 'both');

% Allocate auxiliary variables.
D = zeros(patRad + height + patRad, patRad + width + patRad);
Dint = zeros(1 + size(D, 1), 1 + size(D, 2));
weiAux = zeros(N, 3);
weiNum = N * length(X);
wei = zeros(weiNum, 3);

% Bpad(rows, cols) is Bpad without the winRad pixel frame.
rows = (winRad + 1):(winRad + patRad + height + patRad);
cols = (winRad + 1):(winRad + patRad + width + patRad);

% Compute the linear coordinates of all the pixels in a image.
z = (1:1:N)';

% Initialize the 1st column of "weiAux" with "z".
weiAux(:, 1) = z;

% Compute the (row,col) coordinates of all the pixels in "z".
[y, x] = ind2sub([height, width], z);

% Auxiliary indexes.
rowSelect1 = (1 + patRad + (1 + patRad)):(1 + patRad + (height + patRad));
colSelect1 = (1 + patRad + (1 + patRad)):(1 + patRad + (width + patRad));
rowSelect2 = (1 + patRad + (1 + patRad)):(1 + patRad + (height + patRad));
colSelect2 = (1 + patRad + (1 - patRad - 1)):(1 + patRad + (width - patRad - 1));
rowSelect3 = (1 + patRad + (1 - patRad - 1)):(1 + patRad + (height - patRad - 1));
colSelect3 = (1 + patRad + (1 + patRad)):(1 + patRad + (width + patRad));
rowSelect4 = (1 + patRad + (1 - patRad - 1)):(1 + patRad + (height - patRad - 1));
colSelect4 = (1 + patRad + (1 - patRad - 1)):(1 + patRad + (width - patRad - 1));

% Each row of "wei" will contain, from left to right:
% - the linear coordinate of a 1st pixel,
% - the linear coordinate of a 2nd pixel,
% - the energy of the difference between the two patches centered on the
%   1st and 2nd pixel, respectively.
for k = 1:1:offsetsNum
    
    % We want to compare the patch centered on pixels (y,x) with that at pixel (yp,xp).
    yp = y + Y(k);
    xp = x + X(k);
    
    % Compute the squared difference between Apad and the shifted Bpad.
    D(:, :) = ( Apad - Bpad(rows + Y(k), cols + X(k)) ).^2;
    
    % Compute the integral image of D.
    Dint(:, :) = integralImage(D, 'upright');
    
    % V(y,x) will contain the weight between pixels (y,x) and (yp,xp).
    V = ...
        Dint(rowSelect1, colSelect1) ...
        - Dint(rowSelect2, colSelect2) ...
        - Dint(rowSelect3, colSelect3) ...
        + Dint(rowSelect4, colSelect4);
    
    % Detect those comparisons with pixels (yp,xp) outside image B.
    mask = (yp >= 1) & (yp <= height) & (xp >= 1) & (xp <= width);
    
    % Store the valid pixels (yp,xp) in linear form.
    weiAux(mask, 2) = sub2ind([height, width], yp(mask), xp(mask));
    
    % Compute and store the Gaussian weights.
    weiAux(:, 3) = exp( - V(:) / (intSigma ^2) );
    
    % Remove the non valid weights, i.e. those associated to pixels (yp,xp) outside B.
    weiAux(~mask, 3) = 0;
    
    % Add the weigths assocaited to the offset (Y(k),X(k)).
    wei( (((k - 1) * N) + 1):(k * N), : ) = weiAux;
    
    % Clean weiAux for the next iteration.
    weiAux(:, 2:3) = 0;
    
end

% =========================================================================

if ~QUIET
    
    tElapsed = toc(tStart);
    fprintf('\n>>> nlm <<< time: %f\n', tElapsed);
    
end

end

