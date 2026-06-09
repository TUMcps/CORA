function res = testnn_nnConvTranspose2DLayer_evalutateZonotopeBatch()
% testnn_nnConvTranspose2DLayer_evalutateZonotopeBatch - test 
% nnConvTranspose2DLayer/evalutateZonotopeBatch function
%
% Syntax:
%    res = testnn_nnConvTranspose2DLayer_evalutateZonotopeBatch()
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

% Authors:       Benedikt Kellner
% Written:       01-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Reset random number generator.
rng('default')

% Specify batch size.
bSz = 16;
% Specify input dimensions.
imgSz = [16 16 3];
inDim = prod(imgSz);
% Specify number of filter and size.
kerH = 3;
kerW = 3;
chIn = imgSz(3);
chOut = 10;
% Specify number of generators.
numGen = 10;
% Specify number of random samples for validation.
N = 100;

% Instantiate random layer.
% nnConvTranspose2DLayer expects weights as (H, W, Out, In)
W = rand([kerH kerW chOut chIn]);
b = rand([chOut 1]);
convl = nnConvTranspose2DLayer(W,b);

% Instantiate neural networks with only one layer.
nn = neuralNetwork({convl});
nn.setInputSize(imgSz);
% Prepare the neural network for the batch evaluation.
options = nnHelper.validateNNoptions(struct());
options.nn.train.num_init_gens = numGen;
options.nn.poly_method = 'bounds'; % Needs to be valid for training mode
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
    assert(all(contains(Yi,ysi,'exact',1e-6)));
end

end

% ------------------------------ END OF CODE ------------------------------
