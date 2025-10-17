function res = test_confPredictor_calibrate_05_classOutlier
% test_confPredictor_calibrate_05_classOutlier - unit test for 
%    calibration mit outlier detection of interval and zono-conformal 
%    predictors for classification tasks
%
% Syntax:
%    res = test_confPredictor_calibrate_05_classOutlier
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
% Written:       01-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% user specifications

% number of data points for training/validation, testing and identification
n_cal = 50;

% define predictor model
c1 = [-5; -5]; c2 = [0; -1]; c3 = [3; 1]; % centers for the three classes
f = @(x,u) [sum(10-(x-c1+u(1:2,:)).^2); sum(10-(x-c2).^2+u(3,:)); sum(10-(x-c3).^2+u(4,:))];
dim_x = 2;
dim_u = 4;
dim_y = 3;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create calibration data
data_cal = aux_createData(n_cal, dim_y, [c1, c2, c3]);

% calibration with different outlier detection methods
options.cs.constraints = "gen"; % outlier detection not implemented for halfspace constraints
for type = ["zonotope", "interval"]
    % create predictor
    pred = confPredictor(sys,type,"class");
    for outMethod = ["RMSE", "search", "searchG","MILP"]
        options.cs.outMethod = outMethod;

        % calibrate with default generator template and 1 outlier
        n_out = 1;
        pred = calibrate(pred,data_cal,n_out,options);
        assert(isa(pred.U,type))
        [conservatism1, coverage1] = evaluateData(pred,data_cal);
        if isMATLABReleaseOlderThan('R2024a') && strcmp(outMethod,'MILP') ...
                && ~isfinite(conservatism1)
            % MILP-solver from older Matlab versions does not find feasible
            % solution for this problem
            display("No valid uncertainty set found.")
            continue
        else
            assert(coverage1>=(n_cal-n_out)/n_cal-1e-6);
        end

        % downsize uncertainty set
        U = enlarge(pred.U,0.99);
        pred = setU(pred,U);
        [conservatism1_, coverage1_] = evaluateData(pred,data_cal);
        assert(conservatism1 >= conservatism1_)
        assert(coverage1_ <= coverage1)

        % calibrate with default generator template and 3 outliers
        n_out = 3;
        pred = calibrate(pred,data_cal,n_out,options);
        assert(isa(pred.U,type))
        [conservatism3, coverage3, ~] = evaluateData(pred,data_cal);
        assert(coverage3>=(n_cal-n_out)/n_cal-1e-6);

        % downsize uncertainty set
        U = enlarge(pred.U,0.99);
        pred = setU(pred,U);
        [conservatism3_, coverage3_] = evaluateData(pred,data_cal);
        assert(conservatism3 >= conservatism3_)
        assert(coverage3_ <= coverage3)
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
    data.input(:,i) = 10*rand(2,1)-5+c_i;
    data.target(class_i,i) = 1;
end
end

% ------------------------------ END OF CODE ------------------------------
