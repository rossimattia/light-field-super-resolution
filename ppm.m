
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

function [x, history] = ppm(P, q, r, beta, iters, init)
% PPM solves the following problem using the Proximal Point Method:
%
%   minimize_x     (0.5 * x' * P * x) + (q' * x) + r
%
% INPUT:
% P - a nxn positive semi definite matrix (double).
% q - a nx1 vector (double).
% r - a scalar (double).
% beta - the proximal constant (double).
% iters - the desired number of iterations.
% init - the initial guess for the minimizer (double), or an empty array [].
%
% OUTPUT:
% x - the computed minimizer (a nx1 vector of doubles).
% history - a structure containing the objective value, and the residual
%           norm of x at each iteration.

% =========================================================================

tStart = tic;

% ==== Global constants ===================================================

% QUIET = 1 suppresses any screen output.
QUIET = 1;

% Absolute tolerance.
EPS_ABS = 1e-4;

% ==== PPM solver =========================================================

% Size of P.
n = size(P, 1);

% Allocate the solver variables.
x = zeros(n, 1);
x_old = zeros(n, 1);

% Allocate the struct history.
history.objval = zeros(iters, 1);
history.r_norm = zeros(iters, 1);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\n', 'iter', 'r norm', 'objective');
end

% Compute the threshold for the stopping condition.
threshold = EPS_ABS * sqrt(n);

% Precompute the matrix P + ((1/beta) * I) with I the identity matrix.
A = P + ((1 / beta) * speye(n));

for k = 1:1:iters
    
    % Updates x.
    if (k == 1) && (~isempty(init))
        
        x(:) = init;
    
    else
        
        [x(:), pcgFlag] = pcg(A, (x / beta) - q);
        
        % if (pcgFlag ~= 0)
        %     error('ERROR: pcg did not converge !!!\n\n');
        % end
        
    end
    
    % Check the objective function value.
    history.objval(k)  = objfun(P, q, r, x);
    
    % Compute the residual.
    history.r_norm(k) = norm(x - x_old);

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.2f\n', k, history.r_norm(k), history.objval(k));
    end
    
    % Evaluate the stopping criterion.
    if (history.r_norm(k) < threshold)
        
        break;
    
    end
    
    % Save the current value of x.
    x_old(:) = x;
    
end

if ~QUIET
    tElapsed = toc(tStart);
    fprintf('\n>>> ppm <<< time: %f\n', tElapsed);
end

end

function val = objfun(P, q, r, x)

val = (0.5 * x' * P * x) + (q' * x) + r;

end

