function res = test_confPredictor_evaluateData_01_reg
% test_confPredictor_evaluateData_01_reg - unit test for evaluating set-based 
%   predictors for regression tasks
%
% Syntax:
%    res = test_confPredictor_evaluateData_01_reg
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

for type = ["zonotope", "interval", "conformalAdd"]
    % test zono-conformal and interval predictors
    for constr = ["gen", "half"]
        % test generator and halfspace constraints
        if strcmp(type, "conformalAdd") && strcmp(constr,"half")
            % skip second iteration for conformalAdd
            continue
        end

        % create predictor
        pred = confPredictor(sys,type,"reg");

        % calibrate predictor
        options.cs.constraints = constr;
        pred = calibrate(pred,data_cal,0,options);
        [conservatism, coverage, Y] = evaluateData(pred,data_cal);
        assert(coverage == 1); % coverage over calibration data must be 1
        assert(length(Y) == n_cal)
        if strcmp(type, "conformalAdd")
            assert(isa(Y{1},"interval"))
            assert(isa(Y{end},"interval"))
        else
            assert(isa(Y{1},type))
            assert(isa(Y{end},type))
        end

        % increase size of uncertainty set
        U = enlarge(pred.U,1.01);
        pred1 = setU(pred,U);
        [conservatism1, coverage1, Y1] = evaluateData(pred1,data_cal);
        assert(contains(Y1{1},Y{1}))
        assert(contains(Y1{end},Y{end}))
        assert(coverage1 == 1); % coverage must be 1
        assert(conservatism1 >= conservatism); % conservatism must be bigger

        % decrease size of uncertainty set
        U = enlarge(pred.U,0.99);
        pred2 = setU(pred,U);
        [conservatism2, coverage2, Y2] = evaluateData(pred2,data_cal);
        assert(contains(Y{1},Y2{1}))
        assert(contains(Y{end},Y2{end}))
        assert(coverage2 < 1); % coverage must be smaller than 1
        assert(conservatism2 < conservatism); % conservatism must be smaller

        % compute conservatism and coverage
        vol = 0;
        num_cont = 0;
        for i = 1:n_cal
            vol = vol + volume(Y2{i}); % accumulate volume of each zonotope
            num_cont = num_cont + contains(Y2{i},data_cal.target(:,i));
        end
        assert(withinTol(vol/n_cal,conservatism2,1e-6))
        assert(withinTol(num_cont/n_cal,coverage2,1e-6))
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
