function res = example_confPredictor_calibrate_02_regOutlier
% example_confPredictor_calibrate_02_regOutlier - example for identifying 
%       and calibrating a set-based predictor using different outlier
%       detection methods
%
% Syntax:
%    res = example_confPredictor_calibrate_02_regOutlier
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
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


% number of data points for calibration
n_cal = 100;

% create calibration data
data_cal = aux_createData(n_cal);

% define predictor model 
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(x.*u(1:2,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(x.*u(3:4,:),1)];
dim_x = 2;
dim_u = 4;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create zono-conformal predictor 
pred_zcp = confPredictor(sys,"zonotope","reg");

% create interval-conformal predictor
pred_icp = setType(pred_zcp,"interval");

% create conformal predictor
pred_cp = setType(pred_zcp,"conformalAdd");

% set initial uncertainty sets
G_u = eye(pred_zcp.nrOfUncertainties);
U_init = zonotope(zeros(pred_zcp.nrOfUncertainties, 1),G_u);
pred_zcp = setUinit(pred_zcp, U_init);
pred_icp = setUinit(pred_icp, interval(U_init));

figure; 

% calibration with different outlier detection methods
n_out = 3;

% conformal prediction
pred_cp = calibrate(pred_cp,data_cal,n_out);
[conservatism_cp, coverage_cp, Y_cp] = evaluateData(pred_cp,data_cal);
fprintf("Conservatism CP: %.2f \n", conservatism_cp)
fprintf("Coverage CP: %.2f \n\n", coverage_cp)

% interval and zono-conformal predictors
i = 1;
for outMethod = ["RMSE", "search", "searchG", "MILP"]
    options.cs.outMethod = outMethod;
    pred_zcp = calibrate(pred_zcp,data_cal,n_out,options);
    pred_icp = calibrate(pred_icp,data_cal,n_out,options);

    % evaluate prediction sets
    [conservatism_zcp, coverage_zcp, Y_zcp] = evaluateData(pred_zcp,data_cal);
    [conservatism_icp, coverage_icp, Y_icp] = evaluateData(pred_icp,data_cal);

    % plot predictions for the first data point
    subplot(2,2,i); hold on;
    plot(Y_zcp{1},[1 2], 'DisplayName','ZCP')
    plot(Y_icp{1},[1 2], 'DisplayName','ICP')
    plot(Y_cp{1},[1 2], 'DisplayName','CP')
    plot(data_cal.target(1,1),data_cal.target(2,1),'rx')
    legend;
    title(sprintf('Outlier Detection %s', outMethod))
    fprintf('### Outlier Detection %s ### \n', outMethod)
    fprintf("Conservatism ZCP: %.2f, ICP: %.2f \n", conservatism_zcp, ...
        conservatism_icp)
    fprintf("Coverage ZCP: %.2f, ICP: %.2f \n\n", coverage_zcp, ...
        coverage_icp)
    i = i + 1;
end
res = true;
end


% Auxiliary functions -----------------------------------------------------

function data = aux_createData(n_data)

% define data generation function
f_true = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + x(1,:).*u(1,:); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + x(2,:).*u(2,:)];
n_u = 2;
n_x = 2;
v = randn(n_u,1);

% generate new data
U = 1/5 * zonotope([zeros(n_u,1) ones(n_u,1) 0.1*v]);
X = 5 * interval(-ones(n_x,1),ones(n_x,1));
x = randPoint(X,n_data);
u = randPoint(U,n_data);
y = f_true(x,u);

data.input = x;
data.target = y;
end

% ------------------------------ END OF CODE ------------------------------
