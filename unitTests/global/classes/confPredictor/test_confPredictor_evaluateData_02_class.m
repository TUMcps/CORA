function res = test_confPredictor_evaluateData_02_class
% test_confPredictor_evaluateData_02_class - unit test for evaluating set-based 
%   predictors for classification tasks
%
% Syntax:
%    res = test_confPredictor_evaluateData_02_class
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

% define predictor model 
c1 = [-5; -5]; c2 = [0; -1]; c3 = [3; 1]; % centers for the three classes
f = @(x,u) [sum(10-(x-c1+u(1:2,:)).^2); sum(10-(x-c2).^2+u(3,:)); sum(10-(x-c3).^2+u(4,:))];
dim_x = 2;
dim_u = 4;
dim_y = 3;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% create calibration data
data_cal = aux_createData(n_cal, dim_y, [c1, c2, c3]);

for type = ["zonotope", "interval"]
    % test zono-conformal and interval predictors
    for constr = ["gen", "half"]
        % test generator and halfspace constraints

        % create predictor
        pred = confPredictor(sys,type,"class");

        % calibrate predictor
        options.cs.constraints = constr;
        pred = calibrate(pred,data_cal,0,options);
        [conservatism, coverage, Y] = evaluateData(pred,data_cal);
        assert(coverage == 1); % coverage over calibration data must be 1
        assert(length(Y) == n_cal)
        if strcmp(type,"conformalAdd")
            assert(isnumeric(Y{1}))
            assert(isnumeric(Y{end}))
        else
            assert(isa(Y{1},type))
            assert(isa(Y{end},type))
        end

        % increase size of uncertainty set
        if strcmp(type,"conformalAdd")
            U = pred.U * 1.01;
        else
            U = enlarge(pred.U,1.01);
        end
        pred1 = setU(pred,U);
        [conservatism1, coverage1, Y1] = evaluateData(pred1,data_cal);
        if strcmp(type,"conformalAdd")
            assert(all(contains(string(Y{1}),string(Y1{1}))))
            assert(all(contains(string(Y{end}),string(Y1{end}))))
        else
            assert(contains(Y1{1},Y{1}))
            assert(contains(Y1{end},Y{end}))
        end
        assert(coverage1 == 1); % coverage must be 1
        assert(conservatism1 >= conservatism); % conservatism must be bigger

        % decrease size of uncertainty set
        if strcmp(type,"conformalAdd")
            U = pred.U * 0.99;
        else
            U = enlarge(pred.U,0.99);
        end
        pred2 = setU(pred,U);
        [conservatism2, coverage2, Y2] = evaluateData(pred2,data_cal);
        if strcmp(type,"conformalAdd")
            assert(all(contains(string(Y2{1}),string(Y{1}))))
            assert(all(contains(string(Y2{end}),string(Y{end}))))
        else
            assert(contains(Y{1},Y2{1}))
            assert(contains(Y{end},Y2{end}))
        end
        assert(coverage2 < 1); % coverage must be smaller than 1
        assert(conservatism2 < conservatism); % conservatism must be smaller

        % compute conservatism and coverage
        vol = 0;
        num_correct = 0;
        vec = 1:dim_y;
        for i = 1:n_cal
            class_correct = vec(logical(data_cal.target(:,i)));
            for class = 1:dim_y
                vol = vol + aux_classContains(Y2{i},class); % accumulate number of possible classes
            end
            num_correct = num_correct + aux_classContains(Y2{i},class_correct);
        end
        assert(withinTol(vol/n_cal,conservatism2,1e-6))
        assert(withinTol(num_correct/n_cal,coverage2,1e-6))
    end
end

% test type conformalAdd
% create predictor
pred = confPredictor(sys,"conformalAdd","class");

% calibrate predictor
options.cs.constraints = "half";
pred = calibrate(pred,data_cal,0,options);
[conservatism, coverage, Y] = evaluateData(pred,data_cal);
assert(coverage == 1); % coverage over calibration data must be 1
assert(length(Y) == n_cal)
assert(isnumeric(Y{1}))
assert(isnumeric(Y{end}))

% increase size of uncertainty set
U = pred.U * 1.01;
pred1 = setU(pred,U);
[conservatism1, coverage1, Y1] = evaluateData(pred1,data_cal);
assert(all(contains(string(Y{1}),string(Y1{1}))))
assert(all(contains(string(Y{end}),string(Y1{end}))))

assert(coverage1 == 1); % coverage must be 1
assert(conservatism1 >= conservatism); % conservatism must be bigger

% decrease size of uncertainty set
U = pred.U * 0.99;
pred2 = setU(pred,U);
[conservatism2, coverage2, Y2] = evaluateData(pred2,data_cal);
assert(all(contains(string(Y2{1}),string(Y{1}))))
assert(all(contains(string(Y2{end}),string(Y{end}))))
assert(coverage2 < 1); % coverage must be smaller than 1
assert(conservatism2 < conservatism); % conservatism must be smaller

% compute conservatism and coverage
vol = 0;
num_correct = 0;
vec = 1:dim_y;
for i = 1:n_cal
    class_correct = vec(logical(data_cal.target(:,i)));
    vol = vol + length(Y2{i}); % accumulate number of predicted classes
    num_correct = num_correct + contains(string(class_correct),string(Y2{i}));
end
assert(withinTol(vol/n_cal,conservatism2,1e-6))
assert(withinTol(num_correct/n_cal,coverage2,1e-6))

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
