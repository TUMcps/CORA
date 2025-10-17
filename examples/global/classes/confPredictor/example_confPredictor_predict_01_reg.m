function res = example_confPredictor_predict_01_reg
% example_confPredictor_predict_01_reg - unit test for computing the 
%   predictions of set-based predictors for regression tasks
%
% Syntax:
%    res = example_confPredictor_predict_01_reg
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

% define predictor model
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(x.*u(1:2,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(x.*u(3:4,:),1)];
dim_x = 2;
dim_u = 4;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create and calibrate predictors
pred_zcp = confPredictor(sys,"zonotope","reg");
pred_zcp = calibrate(pred_zcp,data_cal);
pred_icp = confPredictor(sys,"interval","reg");
pred_icp = calibrate(pred_icp,data_cal);
pred_cp = confPredictor(sys,"conformalAdd","reg");
pred_cp = calibrate(pred_cp,data_cal);

% generate and plot prediction sets for the first four data point
figure;
for i = 1:4
    x_i = data_cal.input(:,i);
    y_i = data_cal.target(:,i);
    Y_zcp = predict(pred_zcp,x_i);
    Y_icp = predict(pred_icp,x_i);
    Y_cp = predict(pred_cp,x_i);

    subplot(2,2,i); hold on;
    plot(Y_zcp,[1 2], 'DisplayName','ZCP')
    plot(Y_icp,[1 2], 'DisplayName','ICP')
    plot(Y_cp,[1 2], 'DisplayName','CP')
    plot(y_i(1),y_i(2),'kx')
    xlabel("y_1")
    ylabel("y_2")
    legend;
    title(sprintf('Data point %d', i))
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
