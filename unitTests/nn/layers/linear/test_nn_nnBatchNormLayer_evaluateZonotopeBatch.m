function res = test_nn_nnBatchNormLayer_evaluateZonotopeBatch()
% test_nn_nnBatchNormLayer_evaluateZonotopeBatch - test 
% nnBatchNormLayer/evaluateZonotopeBatch function
%
% Syntax:
%    res = test_nn_nnBatchNormLayer_evaluateZonotopeBatch()
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
% Written:       03-January-2025
% Last update:   21-January-2026 (test moving stats and stored stats)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Reset random number generator.
rng('default')

% Specify batch size.
bSz = 7;
% Specify input and output dimensions.
inSize = [2 3 5];
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Set default options.
options = nnHelper.validateNNoptions(struct());
options.nn.train.num_init_gens = numGen;
options.nn.poly_method = 'bounds';
options.nn.interval_center = true;

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
% Instantiate neural networks with only one layer.
nn = neuralNetwork({batchNorml});
nn.setInputSize(inSize);
% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numGen;
nn.prepareForZonoBatchEval(zeros([prod(inSize) 1]),options);

% Create random batch of input zonotopes.
cx = rand([prod(inSize) 2 bSz]);
Gx = rand([prod(inSize) numGen bSz]);

% Test options that use stored or moving stats. Calculating batch norm 
% stats with sets is not well defined.

% 1. Moving stats ---------------------------------------------------------

% Set option to use the moving stats.
options.nn.batch_norm_stats = 'moving_stats';
% Propagate batch of zonotopes.
[cy,Gy] = nn.evaluateZonotopeBatch(cx,Gx,options);
    
% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output zonotope from the batch.
    Xi = zonotope(cx(:,i),Gx(:,:,i));
    Yi = zonotope(cy(:,i),Gy(:,:,i));
    % Add small interval to avoid numerical errors.
    Yi = Yi + interval(-1e-6*ones([prod(inSize) 1]), ...
        1e-6*ones([prod(inSize) 1]));
    % Sample random points.
    xsi = randPoint(Xi,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi,options);
    % Check if all samples are contained.
    assert(all(contains(Yi,ysi)));
end

% 2. Stored stats ---------------------------------------------------------

% 'Store the stats'.
options.nn.batch_norm_stats = 'calc_stats';
nn.evaluate(xsi,options);

% Set option to use stored stats. 
options.nn.batch_norm_stats = 'stored_stats';
% Propagate batch of zonotopes.
[cy,Gy] = nn.evaluateZonotopeBatch(cx,Gx,options);
    
% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output zonotope from the batch.
    Xi = zonotope(cx(:,i),Gx(:,:,i));
    Yi = zonotope(cy(:,i),Gy(:,:,i));
    % Add small interval to avoid numerical errors.
    Yi = Yi + interval(-1e-6*ones([prod(inSize) 1]), ...
        1e-6*ones([prod(inSize) 1]));
    % Sample random points.
    xsi = randPoint(Xi,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi,options);
    % Check if all samples are contained.
    assert(all(contains(Yi,ysi)));
end

end

% ------------------------------ END OF CODE ------------------------------
