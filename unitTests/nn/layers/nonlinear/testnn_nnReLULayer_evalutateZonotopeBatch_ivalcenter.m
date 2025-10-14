function res = testnn_nnReLULayer_evalutateZonotopeBatch_ivalcenter()
% testnn_nnReLULayer_evalutateZonotopeBatch_ivalcenter - test 
% nnReLULayer/evalutateZonotopeBatch_ivalcenter function with interval center
%
% Syntax:
%    res = testnn_nnReLULayer_evalutateZonotopeBatch_ivalcenter()
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

% Specify evaluation options.
options = nnHelper.validateNNoptions(struct());
options.nn.train.num_init_gens = numGen;
options.nn.poly_method = 'bounds';
options.nn.interval_center = true;

% Instantiate random layer.
relul = nnReLULayer;
% Instantiate neural networks with only one layer.
nn = neuralNetwork({relul});
% Prepare the neural network for the batch evaluation.
nn.prepareForZonoBatchEval(zeros([inDim 1]),options);

% Create random batch of input zonotopes.
cx = sort(rand([inDim 2 bSz]),2);
Gx = rand([inDim numGen bSz]);

% Propagate batch of zonotopes.
[cy,Gy] = nn.evaluateZonotopeBatch(cx,Gx,options);

% Check if all samples are contained.
for i=1:bSz
    % Instantiate i-th input and output zonotope from the batch.
    cxi = 1/2*(cx(:,2,i) + cx(:,1,i));
    dxi = 1/2*(cx(:,2,i) - cx(:,1,i));
    Xi = zonotope(cxi,Gx(:,:,i)) + interval(-dxi,dxi);
    cyi = 1/2*(cy(:,2,i) + cy(:,1,i));
    dyi = 1/2*(cy(:,2,i) - cy(:,1,i));
    Yi = zonotope(cyi,Gy(:,:,i)) + interval(-dyi,dyi);
    % Sample random points.
    xsi = randPoint(Xi,N);
    % Propagate samples.
    ysi = nn.evaluate(xsi);
    % Check if all samples are contained.
    assert(all(contains(Yi,ysi)));
end

end

% ------------------------------ END OF CODE ------------------------------
