function res = test_nn_nnBatchNormLayer_evaluateIntervalBatch()
% test_nn_nnBatchNormLayer_evaluateIntervalBatch - test 
% nnBatchNormLayer/evaluate with a batch of intervals 
%
% Syntax:
%    res = test_nn_nnBatchNormLayer_evaluateIntervalBatch()
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
bSz = 16;
% Specify input and output dimensions.
inSize = [5 7 9];
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Set default options.
options = nnHelper.validateNNoptions(struct());
options.nn.train.num_init_gens = numGen;
options.nn.poly_method = 'bounds';

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
cx = rand([prod(inSize) bSz]);
Gx = rand([prod(inSize) numGen bSz]);
% Compute enclosing intervals.
rx = sum(abs(Gx),2);
xIval = interval(cx - rx(:,:),cx + rx(:,:));

% Test options that use stored or moving stats. Calculating batch norm 
% stats with sets is not well defined.

% 1. Moving stats ---------------------------------------------------------

% Set option to use the moving stats.
options.nn.batch_norm_stats = 'moving_stats';
% Propagate batch of intervals.
yIval = nn.evaluate(xIval,options);
    
% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output interval from the batch.
    xIvali = interval(xIval.inf(:,i),xIval.sup(:,i));
    yIvali = interval(yIval.inf(:,i),yIval.sup(:,i));
    % Sample random points.
    xsi = randPoint(xIvali,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi,options);
    % Check if all samples are contained.
    assert(all(contains(yIvali,ysi)));
end

% 2. Stored stats ---------------------------------------------------------

% 'Store the stats'.
options.nn.batch_norm_stats = 'calc_stats';
nn.evaluate(xsi,options);

% Set option to use stored stats. 
options.nn.batch_norm_stats = 'stored_stats';
% Propagate batch of intervals.
yIval = nn.evaluate(xIval,options);
    
% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output interval from the batch.
    xIvali = interval(xIval.inf(:,i),xIval.sup(:,i));
    yIvali = interval(yIval.inf(:,i),yIval.sup(:,i));
    % Sample random points.
    xsi = randPoint(xIvali,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi,options);
    % Check if all samples are contained.
    assert(all(contains(yIvali,ysi)));
end

end

% ------------------------------ END OF CODE ------------------------------
