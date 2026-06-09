function res = test_nn_nnBatchNormLayer_evaluate()
% test_nn_nnBatchNormLayer_evaluate - test nnBatchNormLayer/evaluate with for 
% moving stats and stored stats.
%
% Syntax:
%    res = test_nn_nnBatchNormLayer_evaluate()
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
% See also: -

% Authors:       Lukas Koller
% Written:       21-January-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Reset random number generator.
rng('default')

% Specify batch size.
bSz = 64;
% Specify input and output dimensions.
inSize = [5 7 9];

% Obtain a random input.
x = rand([inSize bSz]);

% Randomly initialize the stats (mean and variance).
scale = rand([1 1 inSize(3)]);
offset = rand([1 1 inSize(3)]);
movVar = rand(inSize);
movMean = rand(inSize);

% Instantiate random layer.
epsilon = 1e-5;
momentum = 0.99;
batchNorml = nnBatchNormLayer(scale,offset,movVar,movMean, ...
    'batchNorm',epsilon,momentum);

% Initialize a neural network and set the input size.
nn = neuralNetwork({batchNorml});
nn.setInputSize(inSize);

% Set the random stats.
nn.layers{1}.movVar = movVar;
nn.layers{1}.movMean = movMean;

% 1. Moving stats (inference). --------------------------------------------

% Compute the expected output.
x_ = (x - movMean)./sqrt(movVar + epsilon);
y_ = scale.*x_ + offset;

% Compute the default.
y = nn.evaluate(reshape(x,[prod(inSize) bSz]));
% Reshape to a feature map.
y = reshape(y,size(x));

% Check the computed output.
assert(all(withinTol(y,y_,1e-6),'all'));

% Set option to not calculate stats, i.e., use moving stats.
options.nn.batch_norm_stats = 'moving_stats';
% Compute the output.
y = nn.evaluate(reshape(x,[prod(inSize) bSz]),options);
% Reshape to a feature map.
y = reshape(y,size(x));

% Check the computed output.
assert(all(withinTol(y,y_,1e-6),'all'));

% 2. Compute stats (training). --------------------------------------------

% Compute the expected output.
mu = mean(x,[1 2 4]);
sigma = mean((x - mu).^2,[1 2 4]);
x_ = (x - mu)./sqrt(sigma + epsilon);
y_ = scale.*x_ + offset;

% Set option to calculate stats.
options.nn.batch_norm_stats = 'calc_stats';
% Compute the output.
y = nn.evaluate(reshape(x,[prod(inSize) bSz]),options);
% Reshape to a feature map.
y = reshape(y,size(x));

% Check the computed output.
assert(all(withinTol(y,y_,1e-6),'all'));

% 3. Stored stats (robust training). --------------------------------------
% Use the stats that were stored during the previous forward pass.

% Set option to calculate stats.
options.nn.batch_norm_stats = 'stored_stats';
% Compute the output.
y = nn.evaluate(reshape(x,[prod(inSize) bSz]),options);
% Reshape to a feature map.
y = reshape(y,size(x));

% Check the computed output.
assert(all(withinTol(y,y_,1e-6),'all'));

% ------------------------------ END OF CODE ------------------------------
