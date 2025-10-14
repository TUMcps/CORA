function res = test_nn_neuralNetwork_evaluateZonotopeBatch()
% test_nn_neuralNetwork_evaluateZonotopeBatch - unit test function for 
%     neuralNetwork/evaluateZonotopeBatch: sanity check of points
%
% Syntax:
%    res = test_nn_neuralNetwork_evaluateZonotopeBatch()
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
% Written:       29-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Specify number of input and output dimensions.
n0 = 5;
nK = 7;
% Generate a random neural network.
nn = neuralNetwork.generateRandom( ...
    'NrInputs',n0, ...
    'NrOutputs',nK, ...
    'NrLayers',3, ...
    'NrHiddenNeurons',17 ...
);

% Specify a batch size.
bSz = 13;
% Specify number of generators.
numInitGen = 13;
% Specify number of random samples for validation.
N = 100;

% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numInitGen;
numGen = nn.prepareForZonoBatchEval(zeros([n0 1]),options);
% We subtract 10 generators to check if
% nnActivationLayer/evaluateZonotopeBatch correct adds new generators.
numGen = numGen - randi(numGen - numInitGen);

% Create random batch of input zonotopes.
cx = rand([n0 bSz]);
Gx = rand([n0 numGen bSz]);

% Set random batch entries to only be points, i.e., all generators to zero.
zidx = randi(bSz,3);
Gx(:,:,zidx) = 0;

% Propagate batch of zonotopes.
[cy,Gy] = nn.evaluateZonotopeBatch(cx,Gx);

% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output zonotope from the batch.
    Xi = zonotope(cx(:,i),Gx(:,:,i));
    Yi = zonotope(cy(:,i),Gy(:,:,i));
    % Sample random points.
    xsi = randPoint(Xi,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi);
    % Check if all samples are contained.
    assert(all(contains(Yi,ysi)));
end

% ------------------------------ END OF CODE ------------------------------
