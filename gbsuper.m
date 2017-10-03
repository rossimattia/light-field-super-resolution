
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

function ZEst = gbsuper( ...
    ZLR, ZInit, ...
    B, D, ...
    patRad, intSigma, dispMax, ...
    lambda0, lambda1, lambda2, ...
    lb, ub, ...
    innerCycles, outerCycles, ...
    warpMode, ...
    guessFlag)
% GBSUPER super-resolves the input low resolution light field.
%
% INPUT:
% ZLR - the low resolution light field to super-resolve (double [0,1]).
% ZInit - an initial guess on the high resolution light field (double [0,1]).
% B - the blur matrix.
% D - the decimation matrix.
% patRad - the size of the (patRad x patRad) patch used in the graph construction.
% intSigma - the standard deviation value used in the graph construction.
% dispMax - the maximum assumed disparity value.
% lambda0 - the weight of the data fidelity term.
% lambda1 - the weight of the warping term.
% lambda2 - the weight of the graph-based regularizer.
% lb - the pixel value lower bound used in the solver.
% ub - the pixel value upper bound used in the solver.
% innerCycles - the maximum number of iteration to be performed by the solver.
% outerCycles - the number of iterations of the GB super-resolution algorithm.
% warpMode - 0 for the Direct warping matrix construction, and 1 for
%            the SQuare-constraint matrix construction.
% guessFlag - 0 to initialized the solver with a zero light field, and 1 to
%             initialize the solver with ZInit.
%
% OUTPUT:
% ZEst - the super resolved light field (double [0,1]). 

% ==== Solver's fixed parameters ==========================================

beta = 1.0;

% ==== Sub light field dimensions =========================================

% Angular resolution.
vRes = size(ZInit, 1);
hRes = size(ZInit, 2);
M = vRes * hRes;

% Spatial resolution (HR).
yRes = size(ZInit{1, 1}, 1);
xRes = size(ZInit{1, 1}, 2);
N = yRes * xRes;

% ==== Perform super resolution ===========================================

% Precompute the fixed problem terms.
zLR = lf2col(ZLR);

% Start the minimization.
zHist = zeros(N * M, outerCycles);
ZGuess = ZInit;
for k = 1:1:outerCycles
    
    fprintf('ITER. %d\n', k);
    
    % Compute the inter-views graph and the warping matrices.
    fprintf('Computing inter-views graph, and warping ...\n');
    [WA, WB] = intergraph(ZGuess, patRad, intSigma, dispMax);
    
    % Make WA symmetric.
    E = spones(WA);
    mask = and(E, E');
    WA = WA .* mask;
    % Consider switching to "or".
    
    % If the inter-views term is active, then compute L ...
    if (lambda2 > 0)
        
        L = stdlap(WA);
        
    % otherwise, set L to zero.
    else
        
        L = sparse(M * N, M * N);
        
    end
    
    % Compute the objective function terms ...
    
    A = D * B;
    At = A';
    AA = At * A;
    
    sumHAFHAF = sparse(M * N, M * N);
    sumHAFH = sparse(M * N, size(D, 1));
    for u = 1:1:M
        
        colSta = ((u - 1) * N) + 1;
        colEnd = ((u - 1) * N) + N;
        
        % Compute F, the portion of WB (or WA) responsible for the warping
        % of the view Z{u} to all its (NON DIAGONAL) neighboring views.
        switch warpMode
            
            case 'SQ'
                
                F = WB;
                
            case 'DR'
                
                F = WA;
                
            otherwise
                
                error('Invalid warping mode !!!\n\n');
                
        end
        F(:, 1:(colSta - 1)) = 0;
        F(:, (colEnd + 1):end) = 0;
        
        % Normalize F.
        auxDiag = sum(F, 2);
        mask = (auxDiag ~= 0);
        auxDiag(mask) = 1 ./ auxDiag(mask); % warning: division by 0 leads to Inf values.
        F = spdiags(auxDiag, 0, M * N, M * N) * F;
        
        Ft = F';
        
        H = warpmask(F, B, D, vRes, hRes);
        H = blkdiag(H{:});
        HH = H' * H;
        
        sumHAFHAF = sumHAFHAF + (Ft * At * HH * A * F);
        sumHAFH = sumHAFH + (Ft * At * HH);
    
    end
    
    P = 2 * (...
        (lambda0 * AA) ...
        + (lambda1 * sumHAFHAF) ...
        + (lambda2 * L) ...
        );
    
    q = (- 2) * ( ...
        (lambda0 * At) ...
        + (lambda1 * sumHAFH) ...
        ) * zLR;
    
    r = 0;
    
    % If guessFlag == true, initialize the solver with ZInit.
    if guessFlag
        
        zGuess = lf2col(ZGuess);
        
    else
        
        zGuess = [];
        
    end
    
    % Call the solver.
    fprintf('Minimizing the objective function ...\n');
    [aux, ~] = ppm(P, q, r, beta, innerCycles, zGuess);
    
    % Project the estimated light field into [lb,ub]^(N*M).
    aux(aux < lb) = lb;
    aux(aux > ub) = ub;
    
    % Save the estimated light field.
    zHist(:, k) = aux;
    
    % Update the estimate for the next iteration.
    ZGuess = col2lf(zHist(:, k), vRes, hRes, yRes, xRes, 1);

end

% Arrange the light field in a (vRes x hRes) cell array.
ZEst = col2lf(zHist(:, end), vRes, hRes, yRes, xRes, 1);

end

