function res = testnn_nnReLULayer_evalutateZonotopeBatch()
% testnn_nnReLULayer_evalutateZonotopeBatch - test 
% nnReLULayer/evalutateZonotopeBatch function
%
% Syntax:
%    res = testnn_nnReLULayer_evalutateZonotopeBatch()
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
inDim = 2;
outDim = 2;
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Instantiate random layer.
relul = nnReLULayer;
% Instantiate neural networks with only one layer.
nn = neuralNetwork({relul});
% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numGen;
nn.prepareForZonoBatchEval(zeros([inDim 1]),options);

% Create random batch of input zonotopes.
cx = rand([inDim bSz]);
Gx = rand([inDim numGen bSz]);

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

end

% ------------------------------ END OF CODE ------------------------------
