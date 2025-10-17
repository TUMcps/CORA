function res = test_confPredictor_predict_01_reg
% test_confPredictor_predict_01_reg - unit test for computing the 
%   predictions of set-based predictors for regression tasks
%
% Syntax:
%    res = test_confPredictor_predict_01_reg
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
% Written:       01-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% user specifications

% number of data points for training/validation, testing and identification
n_cal = 100;

% create calibration data
data_cal = aux_createData(n_cal);
x_i = data_cal.input(:,1);
y_i = data_cal.target(:,1);

% define predictor model 
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(x.*u(1:2,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(x.*u(3:4,:),1)];
dim_x = 2;
dim_u = 4;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

for type = ["zonotope", "interval", "conformalAdd"]
    % test predictors
    % create predictor
    pred = confPredictor(sys,type,"reg");

    % calibrate predictor
    pred = calibrate(pred,data_cal);
    Y_i = predict(pred,x_i);
    assert(dim_y == dim(Y_i))
    if strcmp(type, "conformalAdd")
        assert(isa(Y_i,"interval"))
    else
        assert(isa(Y_i,type))
    end
    assert(contains(Y_i,y_i))

    % center of the prediction set must be equal to the output of point
    % predictor
    f_pred = f(x_i,zeros(dim_u,1));
    assert(all(withinTol(center(Y_i),f_pred,1e-6)))
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
