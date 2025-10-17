function res = test_confPredictor_predict_02_class
% test_confPredictor_predict_02_class - unit test for computing the 
%   predictions of set-based predictors for classification tasks
%
% Syntax:
%    res = test_confPredictor_predict_02_class
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

% define predictor model 
c1 = [-5; -5]; c2 = [0; -1]; c3 = [3; 1]; % centers for the three classes
f = @(x,u) [sum(10-(x-c1+u(1:2,:)).^2); sum(10-(x-c2).^2+u(3,:)); sum(10-(x-c3).^2+u(4,:))];
dim_x = 2;
dim_u = 4;
dim_y = 3;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create calibration data
data_cal = aux_createData(n_cal, dim_y, [c1, c2, c3]);
x_i = data_cal.input(:,1);
y_i = data_cal.target(:,1);
vec = 1:dim_y;
class_i = vec(logical(y_i));

for type = ["zonotope", "interval", "conformalAdd"]
    % test predictors
    % create predictor
    pred = confPredictor(sys,type,"class");

    % calibrate predictor
    pred = calibrate(pred,data_cal);
    Y_i = predict(pred,x_i);
    if strcmp(type, "conformalAdd")
        assert(isnumeric(Y_i))
        assert(contains(string(class_i),string(Y_i)))
    else
        assert(isa(Y_i,type))
        assert(dim_y == dim(Y_i))
        assert(aux_classContains(Y_i,class_i))
        
        % center of the prediction set must be equal to the output of point 
        % predictor
        f_pred = f(x_i,zeros(dim_u,1));
        assert(all(withinTol(center(Y_i),f_pred,1e-6))) 
    end
end

res = true;
end


% Auxiliary functions -----------------------------------------------------

function data = aux_createData(n_data, dim_y, c)
% sample random points around the center vectors of each class
data.input = zeros(size(c,1),n_data);
data.target = zeros(dim_y,n_data);
for i = 1:n_data
    class_i = randi(dim_y);
    c_i = c(:,class_i);
    data.input(:,i) = 5*rand(2,1)-2.5+c_i;
    data.target(class_i,i) = 1;
end
end

function res = aux_classContains(Y,class_true)
% check if the true class is feasible for a given output set

optionsOptim = optimoptions('linprog','Algorithm', 'interior-point', 'Display','off');
res = false;
if isa(Y,'interval')
    Y = zonotope(Y);
end
eta = size(generators(Y),2);
n_y = dim(Y);
T_i = repmat((1:n_y)==class_true, n_y, 1) - eye(n_y);
problem.f = zeros(1,eta);
problem.Aineq = -T_i*generators(Y);
problem.bineq = T_i*center(Y);
problem.lb = -1*ones(eta,1);
problem.ub = ones(eta,1);
problem.options = optionsOptim;
[~,~,exitflag,~] = CORAlinprog(problem);
if exitflag >0
    res = true;
end
end

% ------------------------------ END OF CODE ------------------------------
