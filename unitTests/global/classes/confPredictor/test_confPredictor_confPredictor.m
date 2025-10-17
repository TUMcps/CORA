function res = test_confPredictor_confPredictor
% test_confPredictor_confPredictor - example for initializing
%       set-based predictors
%
% Syntax:
%    res = test_confPredictor_confPredictor
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

% define predictor model 
regr_func = @(x) [ones(size(x)); x];
f = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + sum(regr_func(x).*u(1:4,:),1); ...
    1./(x(1,:).^2+1)+cos(x(2,:)) + sum(regr_func(x).*u(5:8,:),1)];
dim_x = 2;
dim_u = 8;
dim_y = 2;
sys = nonlinearSysDT(@(x,u) x,1,dim_x,dim_u,f,dim_y);

% % default initialization using function handle
% pred = confPredictor(f); % default
% assert(~isempty(pred.sys))
% assert(pred.type == "zonotope")
% assert(pred.task == "reg")

% default initialization using neural network object
layers = {nnLinearLayer(zeros(64,dim_x),zeros(64,1)); nnReLULayer;
    nnLinearLayer(zeros(64,64),zeros(64,1)); nnReLULayer;
    nnLinearLayer(zeros(dim_y,64),zeros(dim_y,1))};
nn = neuralNetwork(layers);
nn.initWeights();
pred = confPredictor(nn); % default
assert(~isempty(pred.sys))
assert(pred.type == "zonotope")
assert(pred.task == "reg")

% default initialization using nonlinearSysDT object
pred = confPredictor(sys); % default
assert(~isempty(pred.sys))
assert(pred.type == "zonotope")
assert(pred.task == "reg")
pred = pred.setType("interval");
assert(pred.type == "interval")

% create different types of predictors for different tasks
for type = ["zonotope", "interval", "conformalAdd"]
    pred = confPredictor(sys,type);
    assert(pred.type == type)
    assert(pred.task == "reg")
    for task = ["reg", "class"]
        pred = confPredictor(sys,type,task);
        assert(pred.task == task)
        assert(pred.type == type)
    end
end
res = true;
end

% ------------------------------ END OF CODE ------------------------------
