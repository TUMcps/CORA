function res = example_confPredictor_evaluateCoverage
% example_confPredictor_evaluateCoverage - example for computing the
%       coverage, confidence, and n_out for datasets in [1]
%
% Syntax:
%    res = example_confPredictor_evaluateCoverage
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
% Written:       25-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dataset_name = "Energy";

switch dataset_name
    case "Energy"
        % energy dataset
        n_m = 77;
        n_y = 2;
        n_theta = 15;
    case "scm1d"
        % SCM1D dataset
        n_m = 897;
        n_y = 4;
        n_theta = 17;
    case "scm20d"
        % SCM20D
        n_m = 897;
        n_y = 4;
        n_theta = 17;
    case "RF1"
        % RF1
        n_m = 900;
        n_y = 8;
        n_theta = 21;
    case "Photovoltaic"
        % Photovoltaic
        n_m = 864;
        n_y = 4;
        n_theta = 42;
    otherwise
        % other datasets from [1]
        n_m = 1000;
        n_y = 2;
        n_theta = 15;
end

% compute coverage for different numbers of outliers
aux_returnCoverageRow(n_m, n_y, n_theta, dataset_name)
res = true;
end


% Auxiliary functions -----------------------------------------------------

function aux_returnCoverageRow(n_m,n_y, n_u, name)
% coverage for zono-conformal predictor
cov0 = confPredictor.evaluateCoverage([],0.9,n_m,n_u,0);
cov1 = confPredictor.evaluateCoverage([],0.9,n_m,n_u,1);
cov5 = confPredictor.evaluateCoverage([],0.9,n_m,n_u,5);

% coverage for dimension-wise conformal prediction
cov0_CP = confPredictor.evaluateCoverage([],0.9,n_m,n_y,0);
cov1_CP = confPredictor.evaluateCoverage([],0.9,n_m,n_y,1);
cov5_CP = confPredictor.evaluateCoverage([],0.9,n_m,n_y,5);

fprintf("\\emph{%s} & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n", name, ...
    cov0_CP, cov0, cov1_CP, cov1, cov5_CP, cov5);
end

% ------------------------------ END OF CODE ------------------------------
