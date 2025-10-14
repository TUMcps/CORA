function res = testnn_nnBatchNormLayer_evalutateIntervalBatch()
% testnn_nnBatchNormLayer_evalutateIntervalBatch - test 
% nnBatchNormLayer/evaluate with a batch of intervals 
%
% Syntax:
%    res = testnn_nnBatchNormLayer_evalutateIntervalBatch()
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Reset random number generator.
rng('default')

% Specify batch size.
bSz = 16;
% Specify input and output dimensions.
inDim = 5;
outDim = 4;
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Instantiate random layer.
batchNorml = nnBatchNormLayer();
% Instantiate neural networks with only one layer.
nn = neuralNetwork({batchNorml});
% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numGen;
nn.prepareForZonoBatchEval(zeros([inDim 1]),options);

% Create random batch of input zonotopes.
cx = rand([inDim bSz]);
Gx = rand([inDim numGen bSz]);
% Compute enclosing intervals.
rx = sum(abs(Gx),2);
xIval = interval(cx - rx(:,:),cx + rx(:,:));

% Propagate batch of intervals.
yIval = nn.evaluate(xIval);

% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output interval from the batch.
    xIvali = interval(xIval.inf(:,i),xIval.sup(:,i));
    yIvali = interval(yIval.inf(:,i),yIval.sup(:,i));
    % Sample random points.
    xsi = randPoint(xIvali,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi);
    % Check if all samples are contained.
    assert(all(contains(yIvali,ysi)));
end

end

% ------------------------------ END OF CODE ------------------------------
