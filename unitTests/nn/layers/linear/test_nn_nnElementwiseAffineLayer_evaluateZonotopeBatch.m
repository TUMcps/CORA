function res = test_nn_nnElementwiseAffineLayer_evaluateZonotopeBatch()
% test_nn_nnElementwiseAffineLayer_evaluateZonotopeBatch - test 
% nnElementwiseAffineLayer/evaluateZonotopeBatch function
%
% Syntax:
%    res = test_nn_nnElementwiseAffineLayer_evaluateZonotopeBatch()
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
inSize = [2 3 5];
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Set default options.
options = nnHelper.validateNNoptions(struct());
options.nn.train.num_init_gens = numGen;
options.nn.poly_method = 'bounds';

% Randomly initialize the stats (mean and variance).
scale = ones([1 1 inSize(3)]);
offset = zeros([1 1 inSize(3)]);
% Instantiate random layer.
affineLayer = nnElementwiseAffineLayer(scale,offset);
% Instantiate neural networks with only one layer.
nn = neuralNetwork({affineLayer});
nn.setInputSize(inSize);
% Prepare the neural network for the batch evaluation.
options.nn.train.num_init_gens = numGen;
nn.prepareForZonoBatchEval(zeros([prod(inSize) 1]),options);

% Create random batch of input zonotopes.
cx = rand([prod(inSize) 2 bSz]);
Gx = rand([prod(inSize) numGen bSz]);

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
