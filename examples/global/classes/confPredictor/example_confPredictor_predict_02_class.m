function res = example_confPredictor_predict_02_class
% example_confPredictor_predict_02_class - unit test for computing the 
%   predictions of set-based predictors for classification tasks
%
% Syntax:
%    res = example_confPredictor_predict_02_class
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

% create and calibrate predictors
pred_zcp = confPredictor(sys,"zonotope","class");
pred_zcp = calibrate(pred_zcp,data_cal);
pred_icp = confPredictor(sys,"interval","class");
pred_icp = calibrate(pred_icp,data_cal);
pred_cp = confPredictor(sys,"conformalAdd","class");
pred_cp = calibrate(pred_cp,data_cal);

% generate and plot prediction sets for the first four data point
figure;
for i = 1:4
    x_i = data_cal.input(:,i);
    Y_zcp = predict(pred_zcp,x_i);
    Y_icp = predict(pred_icp,x_i);
    Y_cp = predict(pred_cp,x_i);
    ypred_cp = zeros(dim(Y_icp),1);
    ypred_cp(Y_cp) = 1;

    % determine which dimensions to plot
    dim1 = find(data_cal.target(:,i*20)==1);
    if dim1 == size(data_cal.target,1)
        dim2 = dim1-1;
    else 
        dim2 = dim1+1;
    end
    subplot(2,2,i); hold on;
    plot(Y_zcp,[dim1 dim2], 'DisplayName','ZCP')
    plot(Y_icp,[dim1 dim2], 'DisplayName','ICP')
    plot(ypred_cp(dim1), ypred_cp(dim2), 'x', 'DisplayName','CP')
    legend;

    % plot correct classification as green shaded area
    xl = xlim;
    yl = ylim;
    xl = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    yl = [min(yl(1),xl(1)) max(yl(2),xl(2))];
    fill([min(xl(1),yl(1)) max(xl(2),yl(2)) max(xl(2),yl(2))],...
        [min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1))],'green',...
        'FaceAlpha',0.1,'EdgeColor','none', 'DisplayName', 'Correct Classification')
    xlabel("y_" + dim1 + " (correct class)")
    ylabel("y_" + dim2)
    xlim(xl)
    ylim(yl)
    title(sprintf('Data point %d', i))
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
