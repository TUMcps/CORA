function res = example_confPredictor_initializeUncertainties
% example_confPredictor_initializeUncertainties - example for uncertainty
%       placement for set-based predictors
%
% Syntax:
%    res = example_confPredictor_initializeUncertainties
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
% Written:       25-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% user specifications

% number of data points for training/validation, testing and identification
dataset_type ="regArtificial"; % or "regArtificial"
n_id = 100;
n_nn = 1000;

% create synthetic data
dataSpecs = aux_loadDataSpecs(dataset_type);
data_train = aux_createData(dataSpecs, n_nn);
data_cal = aux_createData(dataSpecs, n_id);

% train neural network
net = aux_approxFunction(data_train, dataset_type);

% create zono-conformal predictor 
pred = confPredictor(net,"zonotope","reg");

% select uncertainty placement strategy
figure;
for i_unc = 1:4
    subplot(2,2,i_unc); hold on;
    numUnc = 10*i_unc;
    up_spec.numUnc = numUnc;
    for up_method = ["random","QR","outRandom","outRandom*"]
        if strcmp(up_method,"outRandom*")
            % add random generator to generator matrix G_U
            up_spec.method = "outRandom";
            up_spec.addRandGens = true;
        else
            up_spec.method = up_method;
            up_spec.addRandGens = false;
        end
        pred = initializeUncertainties(pred,up_spec);
        pred = calibrate(pred,data_cal);

        % generate and plot the first prediction set
        Y = predict(pred,data_cal.input(:,1));
        plot(Y,[1 2], 'DisplayName',up_method)
    end
    plot(data_cal.target(1,1),data_cal.target(2,1),'kx')
    xlabel("y_1")
    ylabel("y_2")
    legend;
    title(sprintf('Number uncertainties = %d', numUnc))
end

res = true;
end


% Auxiliary functions -----------------------------------------------------

function dataSpecs = aux_loadDataSpecs(func_type)
% load dataset and dataset parameters

% initialize iuncertainty dimension with 0 (will be set to the true
% uncertainty dimension if the true data generation function is known)
x_add = 0;
u_add = 0;
switch func_type
    case "regArtificial"
        dataSpecs.f_true = @(x,u) [5*sin(x(1,:))+x(2,:).^2 + x(1,:).*u(1,:); ...
            1./(x(1,:).^2+1)+cos(x(2,:)) + x(2,:).*u(2,:)];
        n_y = 2;
        n_u = 2;
        n_x = 2;
        x_lim = 5;
        u_lim = 1/x_lim;
        vseed = 31;

    case "regArtificialSad" %Sadeghi's test function number 3
        dataSpecs.f_true = @(x,u) [3*x(1,:).^3+exp(cos(10*x(2,:)).*cos(5*x(1,:)).^2)+exp(sin(7.5*x(3,:)))+u(1,:); ...
            2*x(1,:).^2+exp(cos(10*x(1,:)).*cos(5*x(2,:)).^2)+exp(sin(7.5*x(3,:).^2))+1.5*u(2,:)];
        n_y = 2;
        n_u = 2;
        n_x = 3;
        x_add = 0.5;
        x_lim = 0.5;
        u_add = 0.5;
        u_lim = 0.5;
        vseed = 30;
end
dataSpecs.n_y = n_y;
dataSpecs.n_u = n_u;
dataSpecs.n_x = n_x;

% synthetic data (will be generated directly before training or
% identification) -> specify number of samples for training and
% validation and input space and uncertainty set
v = randn(RandStream("mt19937ar",'seed', vseed),n_u,1);
dataSpecs.U = u_add+u_lim * zonotope([zeros(n_u,1) ones(n_u,1) 0.1*v]);
dataSpecs.X = x_add+x_lim * interval(-ones(n_x,1),ones(n_x,1));
end

function data = aux_createData(dataSpecs, n_data)
% synthetic dataset: generate new data

x = randPoint(dataSpecs.X,n_data);
u = randPoint(dataSpecs.U,n_data);
y = dataSpecs.f_true(x,u);

data.input = x;
data.target = y;
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
    numNeurons = [size(x_train,1) 64 64 size(y_train,1)];
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

    % Instantiate a neural network.
    nn = neuralNetwork(layers);
    % Initialize the weights of the neural networks.
    nn.initWeights();

    % Specify training options.
    options.nn.train = struct( ...
        'optim',nnAdamOptimizer(0.005,0.9,0.999,1e-8,0,false,[],1,1), ...
        'lr_decay',0.5,'lr_decay_epoch',[5 10 15 20],...
        'max_epoch',20,...
        'mini_batch_size',32,...
        'loss','mse',...
        'shuffle_data','every_epoch' ...
    );
    % Train the neural network.
    nn.train(x_train,y_train,[],[],options);
    % Cast the weights of the neural network back to doubles.
    nn.castWeights(double(1));
end
end


% ------------------------------ END OF CODE ------------------------------
