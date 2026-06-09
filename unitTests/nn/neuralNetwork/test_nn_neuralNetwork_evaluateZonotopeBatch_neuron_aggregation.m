function res = test_nn_neuralNetwork_evaluateZonotopeBatch_neuron_aggregation()
% test_nn_neuralNetwork_evaluateZonotopeBatch_neuron_aggregation - unit 
%   test function for  neuralNetwork/evaluateZonotopeBatch: check neuron
%   aggregation capabilities
%
% Syntax:
%    res = test_nn_neuralNetwork_evaluateZonotopeBatch_neuron_aggregation()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate

% Authors:       Lukas Koller
% Written:       18-February-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Specify number of input and output dimensions.
n0 = 5;
nk = 17;
nK = 7;
K = 3;
% Generate a random neural network.
nn = neuralNetwork.generateRandom( ...
    'NrInputs',n0, ...
    'NrOutputs',nK, ...
    'NrLayers',K, ...
    'NrHiddenNeurons',nk, ...
    'ActivationFun','relu' ...
);

% Specify a batch size.
bSz = 13;
% Specify number of generators.
numInitGen = 13;
% Specify number of random samples for validation.
N = 100;

% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numInitGen;
nn.prepareForZonoBatchEval(zeros([n0 1]),options);

% Create random batch of input zonotopes.
cx = rand([n0 bSz]);
Gx = rand([n0 numInitGen bSz]);

% Specify a neuron aggregation function.
options.nn.neuron_aggregation_fun = @(a,layeri,ci,Gi,co,Go) ...
    aux_numUnstableReLU(aux_numReLU(a,layeri,ci,Gi,co,Go),layeri,ci,Gi,co,Go);
% Propagate batch of zonotopes.
[cy,Gy,a] = nn.evaluateZonotopeBatch(cx,Gx,options);
% Obtain the computed number of ReLU neurons.
numReLU = a.numReLU;
% Obtain the computed number of ReLU neurons.
numUnstable = a.numUnstable;

% Check the size of the aggregation result.
assert(isscalar(numReLU));
assert(numReLU == (K-1)*nk + nK);

% Check the aggregation result.
assert(all(size(numUnstable) == [1 bSz]));
assert(all(numUnstable <= numReLU));


% Generate a random neural network without ReLU.
nn = neuralNetwork.generateRandom( ...
    'NrInputs',n0, ...
    'NrOutputs',nK, ...
    'NrLayers',K, ...
    'NrHiddenNeurons',nk, ...
    'ActivationFun','tanh' ...
);
% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numInitGen;
nn.prepareForZonoBatchEval(zeros([n0 1]),options);

% Specify a neuron aggregation function.
options.nn.neuron_aggregation_fun = @aux_numReLU;
% Propagate batch of zonotopes.
[cy,Gy,a] = nn.evaluateZonotopeBatch(cx,Gx,options);
% Obtain the computed number of ReLU neurons.
numReLU = a.numReLU;

% Check the size of the aggregation result.
assert(isscalar(numReLU));
assert(numReLU == 0);

end


% Auxiliary functions -----------------------------------------------------

function a = aux_numUnstableReLU(a,layeri,ci,Gi,co,Go)
    % Count the number of unstable ReLU neuron in a neural network.

    % Initialize the results field.
    if ~isfield(a,'numUnstable')
        a.numUnstable = 0;
    end

    % Check the type of layer.
    if isa(layeri,'nnReLULayer')
        % Compute the bounds of the current input set.
        r = reshape(sum(abs(Gi),2),size(ci));
        l = ci - r;
        u = ci + r;

        % Count the number of unstable neurons in the current layer.
        numUnstablei = sum(l < 0 & u > 0,1);

        % Aggregate the result in results struct.
        a.numUnstable = a.numUnstable + numUnstablei;
    end
end

function a = aux_numReLU(a,layeri,ci,Gi,co,Go)
    % Count the number of ReLU neuron in a neural network.

    % Initialize the results field.
    if ~isfield(a,'numReLU')
        a.numReLU = 0;
    end

    % Check the type of layer.
    if isa(layeri,'nnReLULayer')
        % Obtain the number of dimensions.
        [nk,~] = size(ci);

        % Aggregate the result in results struct.
        a.numReLU = a.numReLU + nk;
    end
end

% ------------------------------ END OF CODE ------------------------------
