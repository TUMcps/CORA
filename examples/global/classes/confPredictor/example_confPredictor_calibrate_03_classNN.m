function res = example_confPredictor_calibrate_03_classNN
% example_confPredictor_calibrate_03_classNN - example for identifying and
%       calibrating different types of set-based predictors for 
%       classification tasks
%
% Syntax:
%    res = example_confPredictor_calibrate_03_classNN
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

% number of data points for training/validation, testing and identification
dataset_type ="classArtificial1"; % or "classArtificial2"
n_id = 100;
n_nn = 1000;

% create training data
dataSpecs = aux_loadDataSpecs(dataset_type);
% create calibration data
data_cal = aux_createData(dataSpecs, n_id);

% train neural network
data_train = aux_createData(dataSpecs, n_nn);
net = aux_approxFunction(data_train, dataset_type);

% create zono-conformal predictor 
pred_zcp = confPredictor(net,"zonotope","class");

% create interval-conformal predictor
pred_icp = setType(pred_zcp,"interval");

% create conformal predictor 
pred_cp = setType(pred_zcp,"conformalAdd");


% set initial uncertainty sets
G_u = eye(pred_zcp.nrOfUncertainties);
U_init = zonotope(zeros(pred_zcp.nrOfUncertainties, 1),G_u(:,end-10:end));
pred_zcp = setUinit(pred_zcp, U_init);
pred_icp = setUinit(pred_icp, interval(U_init));

% calibration
pred_cp = calibrate(pred_cp,data_cal);
pred_zcp = calibrate(pred_zcp,data_cal);
pred_icp = calibrate(pred_icp,data_cal);

% evaluate prediction sets and plot for the first data point
[conservatism_zcp, coverage_zcp, Y_zcp] = evaluateData(pred_zcp,data_cal);
[conservatism_icp, coverage_icp, Y_icp] = evaluateData(pred_icp,data_cal);
[conservatism_cp, coverage_cp, Y_cp] = evaluateData(pred_cp,data_cal);

% plot prediction sets for the first four data point
figure;
for i = 1:4
    dim1 = find(data_cal.target(:,i*20)==1);
    if dim1 == size(data_cal.target,1)
        dim2 = dim1-1;
    else 
        dim2 = dim1+1;
    end
    subplot(2,2,i); hold on;
    plot(Y_zcp{i*20},[dim1 dim2], 'DisplayName','ZCP')
    plot(Y_icp{i*20},[dim1 dim2], 'DisplayName','ICP')
    legend;
    xl = xlim;
    yl = ylim;
    xl_ = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    yl_ = [min(yl(1),xl_(1)) max(yl(2),xl_(2))];
    fill([min(xl_(1),yl_(1)) max(xl_(2),yl_(2)) max(xl_(2),yl_(2))],...
        [min(xl_(1),yl_(1)) max(xl_(2),yl_(2)) min(xl_(1),yl_(1))],'green',...
        'FaceAlpha',0.1,'EdgeColor','none', 'DisplayName', 'Correct Classification')
    xlabel("y_" + dim1 + " (correct class)")
    ylabel("y_" + dim2)
    xlim(xl)
    ylim(yl)
    title(sprintf('Data point %d', i*20))
end

% print results
fprintf("Conservatism ZCP: %.2f, ICP: %.2f, CP: %.2f \n", conservatism_zcp, ...
    conservatism_icp, conservatism_cp)
fprintf("Coverage ZCP: %.2f, ICP: %.2f, CP: %.2f \n\n", coverage_zcp, ...
    coverage_icp, coverage_cp)

res = true;
end


% Auxiliary functions -----------------------------------------------------

function dataSpecs = aux_loadDataSpecs(func_type)
% load dataset and dataset parameters

% initialize uncertainty dimension with 0 (will be set to the true
% uncertainty dimension if the true data generation function is known)
u_add = 0;
x_lim = 5;
switch func_type
    % choose the task
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
