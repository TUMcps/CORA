function res = test_confPredictor_calibrate_02_regOutlier
% test_confPredictor_calibrate_02_regOutlier - unit test for 
%    calibration mit outlier detection of interval and zono-conformal 
%    predictors for regression tasks
%
% Syntax:
%    res = test_confPredictor_calibrate_02_regOutlier
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
n_cal = 10;

% create calibration data
data_cal = aux_createData(n_cal);

% define predictor model 
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(x.*u(1:2,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(x.*u(3:4,:),1)];
dim_x = 2;
dim_u = 4;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% calibration with different outlier detection methods
for type = ["zonotope", "interval"]
    % create predictor
    pred = confPredictor(sys,type,"reg");
    for outMethod = ["MILP", "RMSE", "search", "searchG"] 
        options.cs.outMethod = outMethod;
        for constr = ["gen", "half"]
            if strcmp(outMethod,"MILP") && strcmp(constr,"half")
                % not implemented yet
                continue
            end
            options.cs.constraints = constr;

            % calibrate with default generator template and 1 outlier
            n_out = 1;
            pred = calibrate(pred,data_cal,n_out,options);
            assert(isa(pred.U,type))
            [conservatism1, coverage1] = evaluateData(pred,data_cal);
            fval1 = pred.results.fval;
            if strcmp(outMethod,"RMSE")
                % detected outliers must not be outside of prediction set
                assert(coverage1>(n_cal-n_out)/n_cal-1e-6);
            else
                % outliers must be outside of prediction sets
                assert(withinTol(coverage1, (n_cal-n_out)/n_cal, 1e-6));
            end

            % downsize uncertainty set
            U = enlarge(pred.U,0.99);
            pred = setU(pred,U);
            [conservatism1_, coverage1_] = evaluateData(pred,data_cal);
            assert(conservatism1 > conservatism1_)
            assert(coverage1_ < coverage1)

            % calibrate with default generator template and 3 outliers
            n_out = 3;
            pred = calibrate(pred,data_cal,n_out,options);
            assert(isa(pred.U,type))

            [conservatism3, coverage3, ~] = evaluateData(pred,data_cal);
            fval3 = pred.results.fval;
            if strcmp(outMethod,"RMSE")
                % detected outliers must not be outside of prediction set
                assert(coverage3>(n_cal-n_out)/n_cal-1e-6);
            else
                % outliers must be outside of prediction sets
                assert(withinTol(coverage3, (n_cal-n_out)/n_cal, 1e-6));
                assert(fval3 < fval1)
            end

            % downsize uncertainty set
            U = enlarge(pred.U,0.99);
            pred = setU(pred,U);
            [conservatism3_, coverage3_] = evaluateData(pred,data_cal);
            assert(conservatism3 > conservatism3_)
            assert(coverage3_ < coverage3)
        end
    end
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
