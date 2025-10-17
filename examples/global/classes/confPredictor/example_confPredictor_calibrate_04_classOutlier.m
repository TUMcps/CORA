function res = example_confPredictor_calibrate_04_classOutlier
% example_confPredictor_calibrate_04_classOutlier - example for identifying 
%       and calibrating a set-based predictor for classification tasks 
%       using different outlier detection methods
%
% Syntax:
%    res = example_confPredictor_calibrate_04_classOutlier
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

% define predictor model
c1 = [-5; -5]; c2 = [0; -1]; c3 = [3; 1]; % centers for the three classes
f = @(x,u) [sum(10-(x-c1+u(1:2,:)).^2); sum(10-(x-c2).^2+u(3,:)); sum(10-(x-c3).^2+u(4,:))];
dim_x = 2;
dim_u = 4;
dim_y = 3;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create calibration data
data_cal = aux_createData(n_cal, dim_y, [c1, c2, c3]);

% create zono-conformal predictor 
pred_zcp = confPredictor(sys,"zonotope","class");

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
    dim1 = find(data_cal.target(:,i*20)==1);
    if dim1 == size(data_cal.target,1)
        dim2 = dim1-1;
    else 
        dim2 = dim1+1;
    end
    plot(Y_zcp{i*20},[dim1 dim2], 'DisplayName','ZCP')
    plot(Y_icp{i*20},[dim1 dim2], 'DisplayName','ICP')
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

% ------------------------------ END OF CODE ------------------------------
