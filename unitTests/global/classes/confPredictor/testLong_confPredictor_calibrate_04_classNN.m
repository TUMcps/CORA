function res = testLong_confPredictor_calibrate_04_classNN
% testLong_confPredictor_calibrate_04_classNN - unit test for calibration of
%    set-based predictors for classification tasks
%
% Syntax:
%    res = testLong_confPredictor_calibrate_04_classNN
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

%% user specifications

% number of data points for training/validation, testing and identification
dataset_type ="classArtificial1"; % or "classArtificial2"
n_cal = 100;
n_nn = 1000;

% create synthetic data
dataSpecs = aux_loadDataSpecs(dataset_type);
data_train = aux_createData(dataSpecs, n_nn);
data_cal = aux_createData(dataSpecs, n_cal);

% train neural network
net = aux_approxFunction(data_train, dataset_type);

% create zono-conformal predictor 
pred = confPredictor(net,"zonotope","class");
% use default generator template 
pred1 = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred1,data_cal);
assert(coverage == 1); % coverage over calibration data must be 1
% specify calibration options
options.cs.constraints = 'half';
pred2 = calibrate(pred,data_cal,0,options);
[~, coverage, ~] = evaluateData(pred2,data_cal);
assert(coverage == 1);
assert(withinTol(pred1.results.fval,pred2.results.fval,1e-6))
% set initial uncertainty set
dim_u = pred.nrOfUncertainties;
U_init = zonotope(zeros(dim_u,1), randn(dim_u));
pred = setUinit(pred,U_init);
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isa(pred.U,'zonotope'))

% create interval predictor model
pred = confPredictor(net,"interval","class");
% use default generator template 
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isa(pred.U,'interval'))

% create conformal predictor
pred = confPredictor(net,"conformalAdd","class");
% use default generator template 
pred = calibrate(pred,data_cal);
[~, coverage, ~] = evaluateData(pred,data_cal);
assert(coverage == 1);
assert(isnumeric(pred.U))

res = true;
end


% Auxiliary functions -----------------------------------------------------

function dataSpecs = aux_loadDataSpecs(func_type)
% load dataset and dataset parameters

% initialize iuncertainty dimension with 0 (will be set to the true
% uncertainty dimension if the true data generation function is known)
u_add = 0;
x_lim = 5;
switch func_type
    case "classArtificial1"
        % classification task 1
        dataSpecs.f_true = cell(3,1);
        dataSpecs.f_true{1} = @(x,u) [3*sin(x(1,:)) + u(1,:)];
        dataSpecs.f_true{2} = @(x,u) [x(1,:).^2 + u(1,:)];
        dataSpecs.f_true{3} = @(x,u) [2*x(1,:)-10 + u(1,:)];
        n_y = 3;
        n_u = 1;
        n_x = 2;
        u_lim = 2;

    case "classArtificial2"
        % classification task 2
        dataSpecs.f_true = cell(4,1);
        dataSpecs.f_true{1} = @(x,u) [x(2,:).*sin(x(1,:))+u(1,:)];
        dataSpecs.f_true{2} = @(x,u) [x(1,:).^2+x(2,:)+2*u(1,:)];
        dataSpecs.f_true{3} = @(x,u) [2*x(1,:)-10+x(1,:).*x(2,:)+0.5*u(1,:).^2];
        dataSpecs.f_true{4} = @(x,u) [2*x(1,:)-16+x(2,:).*u(1,:)];
        n_y = 4;
        n_x = 3;
        n_u = 1;
        u_lim = 1;
end
dataSpecs.n_y = n_y;
dataSpecs.n_u = n_u;
dataSpecs.n_x = n_x;

% synthetic data (will be generated directly before training or
% identification) -> specify number of samples for training and
% validation and input space and uncertainty set
v = randn(RandStream("mt19937ar",'seed', 1),n_u,1);
dataSpecs.U = u_add+u_lim * zonotope([zeros(n_u,1) ones(n_u,1) 0.1*v]);
dataSpecs.X = x_lim * interval(-ones(n_x-1,1),ones(n_x-1,1));
end

function data = aux_createData(dataSpecs, n_data)
% synthetic dataset: generate new data

xy = [];
target = [];
for i = 1:length(dataSpecs.f_true)
    % create data points from each class
    n = ceil(n_data/length(dataSpecs.f_true));
    x = randPoint(dataSpecs.X,n);
    u = randPoint(dataSpecs.U,n);
    xy = [xy [x; dataSpecs.f_true{i}(x,u)]];
    target_i = 1:length(dataSpecs.f_true) == i;
    target = [target repmat(target_i', 1,n)];
end
data.input = xy;
data.target = target;
end

function nn = aux_approxFunction(data_nnTrain, dataset_type)
% train neural network

nn_name = sprintf("nnConverted_%s.onnx", dataset_type);
x_train = data_nnTrain.input;
y_train = data_nnTrain.target;

%check if network already exists
try
    nn = neuralNetwork.readONNXNetwork(nn_name);
catch ME
    if ~(strcmp(ME.identifier,'nnet_cnn_onnx:onnx:FileNotFound'))
        rethrow(ME)
    end

    % Specify the number of neurons.
    numNeurons = [size(x_train,1) 32 32 size(y_train,1)];
    % Create the layers of the neural network.
    layers = {};
    for i=1:length(numNeurons)-2
        % Obtain the number of input and output neurons.
        nin = numNeurons(i);
        nout = numNeurons(i+1);
        % Append a linear layer.
        layers{end+1} = nnLinearLayer(zeros(nout,nin),zeros(nout,1));
        % Append a tanh layer.
        layers{end+1} = nnTanhLayer;
    end
    nin = numNeurons(end-1);
    nout = numNeurons(end);
    % Append an output linear layer.
    layers{end+1} = nnLinearLayer(zeros(nout,nin),zeros(nout,1));
    % Append a softmax layer.
    layers{end+1} = nnSoftmaxLayer;

    % Instantiate a neural network.
    nn = neuralNetwork(layers);
    % Initialize the weights of the neural networks.
    nn.initWeights();

    % Specify training options.
    options.nn.train = struct( ...
        'optim',nnAdamOptimizer(0.01,0.9,0.999,1e-8,0,false,[],1,1), ...
        'lr_decay',0.5,'lr_decay_epoch',[5 10 15 20],...
        'max_epoch',20,...
        'mini_batch_size',32,...
        'loss','softmax+log',...
        'shuffle_data','every_epoch' ...
    );
    % Train the neural network.
    nn.train(x_train,y_train,[],[],options);
    % Cast the weights of the neural network back to doubles.
    nn.castWeights(double(1));
end
end

% ------------------------------ END OF CODE ------------------------------
