
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

function [WA, WB] = intergraph(V, patRad, intSigma, dispMax)
% INTERGRAPH computes the graph adjacency matrix used in the graph-based
% regularizer and a graph adjacency matrix containing the square-constraint
% warping matrices as sub matrices.
%
% INPUT:
% V - a light field.
% patRad - the size of the (patRad x patRad) patch to be used in both the
%          graph and warping matrices construction.
% intSigma - the standard deviation value used in the graph construction.
% dispMax - the maximum assumed disparity value.
%
% OUTPUT:
% WA - the graph adjacency matrix associated to the regularizer.
% WB - the graph adjacency matrix associated to the square-constraint warping matrices.
%
% Neither WA nor WB is made symmetric.

% ==== Light field parameters =============================================

% Angular resolution.
vRes = size(V, 1);
hRes = size(V, 2);
M = vRes * hRes;

% Spatial resolution.
yRes = size(V{1, 1}, 1);
xRes = size(V{1, 1}, 2);
N = yRes * xRes;

% ==== Compute the weights ================================================

% We process the views in V sequentially. We center a 3x3 window on the
% central view, and we compute the weights between the central view and the
% other 8 views in the window (the neighboring views).

% Compute the coordinates for the views in the window, with the current
% view (the one that we are processing) in the middle, at (0,0).
[S, T] = meshgrid((-1):1:1);
S = S(:);
T = T(:);

% Remove the coordinate (0,0).
mask = ((S == 0) & (T == 0));
S = S(~mask);
T = T(~mask);

% In the graph WA, each pixel in a view will have at most 2 connections to each
% one of the 4 horizontal/vertical neighboring views, and at most 4 connections
% to each one of the 4 diagonal neighboring views. At most 24 connections in total.
weiA = zeros(N * 24 * M, 3);
counterA = 0;

% In graph WB, each pixel in a view will have at most 2 adjacent connections
% to the 4 horizontal/vertical neighboring views, and the connections lie on
% a (virtual) square, ideally defined by disparity. At most 8 connections in total.
weiB = zeros(N * 8 * M, 3);
counterB = 0;

% Note that WB (like WA) is tiled with NxN matrices WBij that record the
% connections between a pixel in view i (the rows of WBij) and a pixel in
% view j (the columns of Wij). Normalizing the rows of matrix WBij, makes
% it a warping matrix from view j to view i.

% The total number of integer disparities. 
dispNum = (2 * dispMax) + 1;

% dispScore(i,j) will contain the sum of the weights from the i-th pixel
% of V{t,s} to the pixels of the neighboring views in the square associated
% to disparity j.
dispScore = zeros(N, dispNum);

% dispScoreCounter(i,j) stores the number of weights that contributed to
% dispScore(i,j), for normalization purposes.
dispScoreCounter = zeros(N, dispNum);

% Auxiliary linear array.
auxLin = (1:1:N)';

% Start computing WA and WB ...
for u = 1:1:M
    
    % The number of neighboring views.
    neigNum = length(S);
    
    % Initialize the cell array for the weights storage.
    % >>> Consider initializing each cell with a fixed size vector, and a counter !!!
    neigWei = cell(neigNum, 1);
    
    % Compute the coordinates of the view we are processing.
    [t, s] = ind2sub([vRes, hRes], u);
    
    % Compute WA connections from V{t,s} to its neighborhood, and the
    % square-constraint disparity ranges for the pixels in V{t,s}.
    for k = 1:1:neigNum
        
        % We want to compare views V{t,s} and V{tp,sp}.
        tp = t + T(k);
        sp = s + S(k);
        
        % Note that (tp,sp) cannot be (t,s), cause the pair (T(k),S(k))
        % equal to (0,0) has been removed.
        
        % If (tp,sp) is a valid coordinate, then ...
        if (tp >= 1) && (tp <= vRes) && (sp >= 1) && (sp <= hRes)
            
            % Compute the linear coordinate of (tp,sp).
            up = sub2ind([vRes, hRes], tp, sp);
            
            % If V{tp,sp} is NOT a DIAGONAL neighboring view, then ...
            if ((tp == t) || (sp == s))
                
                % Compute the searching window: a (1 x dispNum) window
                % for horizontal neighboring views, and (dispNum x 1) for
                % vertical ones.
                offset.X = (- sign(S(k))) * ((- dispMax):1:dispMax)';
                offset.Y = (- sign(T(k))) * ((- dispMax):1:dispMax)';
                % Note that a pixel of V{t,s} has a disparity d with reference
                % to the left view, but -d with reference to the right one.
                % Multiplication by -sign(S(k)) and -sign(T(k)) is made
                % necessary by this fact.
                
                % In WA each pixel has (at most) 2 connections to each
                % horizontal/vertical neighboring views.
                cntNum = 2;
                
            % otherwise ...    
            else
                
                % Create a mask to select diagonals -1, 0, and 1.
                % mask = triu(ones(dispNum), -1) - triu(ones(dispNum), 2);
                % mask = logical(mask);
                mask = (ones(dispNum) ~= 0);
                
                % Compute the searching window: a diagonal window with the
                % shape defined by "mask".
                [X, Y] = meshgrid((- dispMax):1:dispMax);
                X = X(:);
                Y = Y(:);
                offset.X = (- sign(S(k))) * X(mask(:));
                offset.Y = (- sign(T(k))) * Y(mask(:));
                
                % In WA each pixel has 4 connections to each diagonal
                % neighboring view.
                cntNum = 4;
                
            end
            
            % Compute the weights between the two views.
            neigWei{k} = nlm(V{t, s}, V{tp, sp}, patRad, intSigma, offset);
            
            % Extract the cntNum best matches for each pixel in V{t,s} ...
            auxMax = reshape(neigWei{k}(:, 3), N, []);
            for z = 1:1:cntNum
                
                % Find the best (highest weight) connection for each pixel.
                [~, idx] = max(auxMax, [], 2);
                
                % Set to zero the found connections, in order to make them
                % no longer valid for selection at the next iteration.
                auxMax(sub2ind(size(auxMax), auxLin, idx)) = 0;
                
                % Extract the best connection for each pixel.
                auxWei = ...
                    neigWei{k}( ...
                    ((idx - 1) * N) + auxLin, ...
                    : ...
                    );
                
                % Detect the non valid connections. These may be connections
                % to pixels outside the view support, or connections already
                % selected in a previous iteration.
                mask = (auxWei(:, 3) ~= 0);
                maskNum = sum(mask);
                
                % Store the best (and valid) connections in weiA.
                weiA((counterA + 1):(counterA + maskNum), :) = ...
                    cat(2, ...
                    ((u - 1) * N) + auxWei(mask, 1), ...
                    ((up - 1) * N) + auxWei(mask, 2), ...
                    auxWei(mask, 3) ...
                    );
                counterA = counterA + maskNum;
                
            end
            
            % If V{tp,sp} is NOT a DIAGONAL neighboring view, then ...
            if ((tp == t) || (sp == s))
                
                % Update the square-constraint scores ...
                auxWei = reshape(neigWei{k}(:, 3), N, []);
                mask = (auxWei ~= 0);
                dispScore(mask) = dispScore(mask) + auxWei(mask);
                dispScoreCounter(mask) = dispScoreCounter(mask) + 1;
                
            end
            
        end
        
    end
    
    % Normalize the square-constraint scores.
    mask = (dispScoreCounter ~= 0);
    dispScore(mask) = dispScore(mask) ./ dispScoreCounter(mask);
    
    % Compute the square-constraint disparity range at each pixel.
    [idx1, idx2] = sqconstr(dispScore);
    
    % Resets dispScore and dispScoreCounter, for the view V{t,s} to
    % process at the next iteration.
    dispScore = dispScore * 0;
    dispScoreCounter = dispScoreCounter * 0;
    
    % Compute the warping matrix from V{t,s} to each view in its neighborhood ...
    for k = 1:1:neigNum
        
        % We want to build the warping matrix from V{tp,sp} to V{t,s}.
        tp = t + T(k);
        sp = s + S(k);
        
        % If view (tp,sp) is valid (hence neigWei{k} is not empty) and it
        % is NOT a DIAGONAL view, then ...
        if (~isempty(neigWei{k})) && ((tp == t) || (sp == s))
            
            % Compute the linear coordinate of (tp,sp).
            up = sub2ind([vRes, hRes], tp, sp);
            
            % Extract the weights associated to idx1 for each pixel.
            aux1 = neigWei{k}( ...
                ((idx1 - 1) * N) + auxLin, ...
                : ...
                );
            
            % Extract the weights associated to idx2 for each pixel.
            aux2 = neigWei{k}( ...
                ((idx2 - 1) * N) + auxLin, ...
                : ...
                );
            
            % Normalize the (at most) 2 weights associated to each pixel.
            % auxNorm = (aux1(:, 3) + aux2(:, 3));
            % mask = (auxNorm ~= 0);
            % aux1(mask, 3) = aux1(mask, 3) ./ auxNorm(mask);
            % aux2(mask, 3) = aux2(mask, 3) ./ auxNorm(mask);
            
            % Detect the non zero weights.
            mask1 = (aux1(:, 3) ~= 0);
            mask2 = (aux2(:, 3) ~= 0);
            
            % Compute the number of nonzero weights for each pixel.
            mask1Num = sum(mask1);
            mask2Num = sum(mask2);
            
            % Store the non zero weigths associated to idx1.
            weiB((counterB + 1):(counterB + mask1Num), :) = ...
                cat(2, ...
                ((u - 1) * N) + aux1(mask1, 1), ...
                ((up - 1) * N) + aux1(mask1, 2), ...
                aux1(mask1, 3) ...
                );
            counterB = counterB + mask1Num;
            
            % Store the non-zero weights associated to idx2.
            weiB((counterB + 1):(counterB + mask2Num), :) = ...
                cat(2, ...
                ((u - 1) * N) + aux2(mask2, 1), ...
                ((up - 1) * N) + aux2(mask2, 2), ...
                aux2(mask2, 3) ...
                );
            counterB = counterB + mask2Num;
            
        end
        
    end

end

% Remove empty cells.
weiA(counterA:end, :) = [];
weiB(counterB:end, :) = [];

% Build WA.
WA = sparse( ...
    weiA(:, 1), ...
    weiA(:, 2), ...
    weiA(:, 3), ...
    M * N, M * N ...
    );

% Build WB.
WB = sparse( ...
    weiB(:, 1), ...
    weiB(:, 2), ...
    weiB(:, 3), ...
    M * N, M * N ...
    );

end


function [idx1, idx2] = sqconstr(dispScore)

dispNum = size(dispScore, 2);
% dispMax = (dispNum - 1) / 2;

[~, idx1] = max(dispScore, [], 2);
% disp1 = idx1 - (dispMax + 1);

idxLeft = max(1, idx1 - 1);
idxRight = min(dispNum, idx1 + 1);

maskLeft = (idxLeft == idx1);
maskRight = (idxRight == idx1);

auxLin = (1:1:length(idx1))';
left = sub2ind(size(dispScore), auxLin, idxLeft);
right = sub2ind(size(dispScore), auxLin, idxRight);

diff = dispScore(left) - dispScore(right);
mask = (diff > 0);

idx2 = zeros(size(idx1));
idx2(mask) = idxLeft(mask);
idx2(~mask) = idxRight(~mask);
idx2(maskLeft) = idxRight(maskLeft);
idx2(maskRight) = idxLeft(maskRight);

% disp2 = zeros(size(disp1));
% disp2(mask) = idxLeft(mask) - (dispMax + 1);
% disp2(~mask) = idxRight(mask) - (dispMax + 1);
% disp2(maskLeft) = idxRight(maskLeft) - (dispMax + 1);
% disp2(maskRight) = idxLeft(maskRight) - (dispMax + 1);

end

