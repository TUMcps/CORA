function res = test_confPredictor_evaluateCoverage
% test_confPredictor_evaluateCoverage - unittest for computing the
%       coverage, confidence, and outlier budget
%
% Syntax:
%    res = test_confPredictor_evaluateCoverage
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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

% Authors:       Laura Luetzow
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i = 1:5
    % test with random values
    n_m = randi([800 1200]);
    n_theta = randi([1 10]);
    coverage = 0.8 + 0.15*rand(1);
    confidence = 0.5 + 0.4*rand(1);

    % compute outlier budget
    n_out = confPredictor.evaluateCoverage(coverage, confidence, n_m, n_theta, []);
    if n_out < 0
        % coverage and confidence levels are not reachable for the
        % predictor
        continue
    end

    % compute confidence for the estimated outlier budget and check if
    % higher than target confidence
    confidence_ = confPredictor.evaluateCoverage(coverage, [], n_m, n_theta, n_out);
    assert(confidence_ >= confidence)

    % compute confidence for a bigger outlier budget and check if
    % smaller than target confidence
    confidence_ = confPredictor.evaluateCoverage(coverage, [], n_m, n_theta, n_out+1);
    assert(confidence_ < confidence)

    % compute coverage for the estimated outlier budget and check if
    % higher than target coverage
    coverage_ = confPredictor.evaluateCoverage([], confidence, n_m, n_theta, n_out);
    assert(coverage_ >= coverage)

    % compute coverage for a bigger outlier budget and check if
    % smaller than target coverage
    coverage_ = confPredictor.evaluateCoverage([], confidence, n_m, n_theta, n_out+1);
    assert(coverage_ < coverage)
end

res = true;
end


% ------------------------------ END OF CODE ------------------------------
