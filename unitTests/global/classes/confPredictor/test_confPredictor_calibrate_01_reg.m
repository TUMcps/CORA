function res = test_confPredictor_calibrate_01_reg
% test_confPredictor_calibrate_01_reg - unit test for calibration of
%    set-based predictors for regression tasks
%
% Syntax:
%    res = test_confPredictor_calibrate_01_reg
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
% See also: setPredictor

% Authors:       Laura Luetzow
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% user specifications

% number of data points for training/validation, testing and identification
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
pred = confPredictor(sys,"zonotope","reg");
% use default generator template 
pred1 = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred1,data_cal);
assert(coverage == 1); % coverage over calibration data must be 1
% specify calibration options
options.cs.constraints = 'half';
pred2 = calibrate(pred,data_cal,0,options);
[~, coverage, ~] = evaluateData(pred2,data_cal);
assert(coverage == 1);
assert(all(withinTol(pred1.U.Z,pred2.U.Z,1e-6),'all'))
% set initial uncertainty set
U_init = zonotope(zeros(dim_u,1), randn(dim_u));
pred = setUinit(pred,U_init);
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isa(pred.U,'zonotope'))

% create interval predictor model
pred = confPredictor(sys,"interval","reg");
% use default generator template 
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isa(pred.U,'interval'))

% create conformal predictor
pred = confPredictor(sys,"conformalAdd","reg");
% use default generator template 
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isa(pred.U,'interval'))

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
