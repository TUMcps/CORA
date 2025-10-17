function output = evaluateCoverage(coverage, confidence, n_m, n_theta, n_out)
% evaluateCoverage - Computes one of three parameters: coverage, 
%       confidence, or n_out
%
% Syntax:
%   [coverage] = confPredictor.evaluateCoverage([], confidence, n_m, n_theta, n_out)
%   [confidence] = confPredictor.evaluateCoverage(coverage, [], n_m, n_theta, n_out)
%   [n_out] = confPredictor.evaluateCoverage(coverage, confidence, n_m, n_theta, [])
%
% Inputs:
%   coverage - target probability 1-epslion, such that probability 
%       P(y\inY(x)) \geq 1-epsilon
%   confidence - confidence level 1-zeta, such that probability
%       P(P(y\inY(x)) \geq 1-epsilon) \geq 1-zeta
%   n_m - number of data points
%   n_theta - number of parameters
%   n_out - number of removed data points
%
% Outputs:
%   output - the value of the parameter that was not provided
%
% References:
%    [1] Campi et al. 2009, "Interval predictor models: Identification and 
%        reliability".
%    [2] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based 
%        Uncertainty Quantification for Regression and Classification 
%        Tasks," arXiv, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: confPredictor

% Authors:       Michael Eichelbeck, Laura Luetzow
% Written:       25-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Check which parameter is missing
if isempty(coverage)
    % Solve for coverage
    output = aux_solve_for_coverage(confidence, n_m, n_theta, n_out);
elseif isempty(confidence)
    % Solve for confidence
    output = aux_solve_for_confidence(coverage, n_m, n_theta, n_out);
elseif isempty(n_out)
    % Solve for n_out
    output = aux_solve_for_n_out(coverage, confidence, n_m, n_theta);
else
    throw(CORAerror("CORA:specialError",'Exactly one of coverage, confidence, or n_out must be empty.'));
end
end


% Auxiliary functions -----------------------------------------------------

function coverage = aux_solve_for_coverage(confidence, n_m, n_theta, n_out)
% Uses a numerical approach to find coverage

    function conf = confidence_from_coverage(rel)
        conf = aux_solve_for_confidence(rel, n_m, n_theta, n_out);
    end

    % Binary search to find coverage
    cov_low = 0.0;  % Lower bound
    cov_high = 1.0; % Upper bound
    
    max_iter = 500;
    tol = 1e-6;
    
    for i = 1:max_iter
        rel_mid = (cov_low + cov_high) / 2;
        conf_mid = confidence_from_coverage(rel_mid);
        
        if abs(conf_mid - confidence) < tol
            coverage = rel_mid;
            return;
        elseif conf_mid < confidence
            cov_high = rel_mid;
        else
            cov_low = rel_mid;
        end
    end
    coverage = (cov_low + cov_high) / 2;
end

function confidence = aux_solve_for_confidence(eta, n_m, n_theta, n_out)
    epsilon = 1 - eta;
    %n_m = vpa(n_m);
    
    zeta = 0;
    for i=0:n_out+n_theta-1
        zeta = zeta + nchoosek(n_m,i) * epsilon^i * (1-epsilon)^(n_m-i);
    end
    zeta = nchoosek(n_out+n_theta-1,n_out) * zeta;
    
    confidence = double(1 - zeta);
end

function n_out = aux_solve_for_n_out(coverage, confidence, n_m, n_theta)
    % Find largest k that still gives confidence >= target
    k = 0;
    while true
        current_conf = aux_solve_for_confidence(coverage, n_m, n_theta, k);
        if current_conf < confidence
            n_out = k-1;
            return;
        end
        k = k + 1;
        
        % Safety check to prevent infinite loop
        if k > n_m
            throw(CORAerror('CORA:specialError','Could not find suitable n_out within reasonable range'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
