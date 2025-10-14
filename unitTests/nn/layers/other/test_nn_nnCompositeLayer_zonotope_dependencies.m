function [res] = test_nn_nnCompositeLayer_zonotope_dependencies()
% test_nn_nnCompositeLayer_zonotope_dependencies - test the dependencies
%   preserved by a zonotope propagation through a nnCompositeLayer.
%
% Syntax:
%    res = test_nn_nnCompositeLayer_zonotope_dependencies()
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
% See also: none

% Authors:       Lukas Koller
% Written:       04-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Specify the number of input dimensions.
n0 = 10;
% Specify number of initial generators.
q0 = 7;

% Instantiate layers: one is identity and the other flips the sign.
layers = {
    {nnLinearLayer(eye(n0),0)};
    {nnLinearLayer(-eye(n0),0)}
};

% 1. ADD
% Instantiate a composite layer which adds the computation paths.
addLayer = nnCompositeLayer(layers,'add');
% Instantiate a sum layer that sum all dimensions.
addSumLayer = nnLinearLayer(ones([1 n0]),0);
% Instantiate a neural network.
nnAdd = neuralNetwork({addLayer, addSumLayer});
nnAdd.setInputSize([n0 1]);

% The expected output is 0, because the computation paths cancel each
% other. The cancelation is only observed for point and zonotope
% propagation, not for intervals.
% Check intervals.
Ix = interval(-ones([n0 1]),ones([n0 1]));
Iy = nnAdd.evaluate(Ix);
assert(Iy.inf == -2*n0 & Iy.sup == 2*n0);
% Check points.
xs = randPoint(Ix,100);
ys = nnAdd.evaluate(xs);
assert(all(ys == 0,'all'));
% Check zonotopes.
options.nn.train.num_init_gens = n0;
nnAdd.prepareForZonoBatchEval(xs,options);
% Specify a batch size.
bSz = 3;
[cy,Gy] = nnAdd.evaluateZonotopeBatch( ...
    zeros([n0 bSz]),repmat(eye(n0,q0),1,1,bSz));
assert(all(cy == 0,'all') && all(Gy == 0,'all'));

% 2. ADD
% Instantiate a composite layer which adds the computation paths.
concatLayer = nnCompositeLayer(layers,'concat');
% Instantiate a sum layer that sum all dimensions.
concatSumLayer = nnLinearLayer(ones([1 2*n0]),0);
% Instantiate a neural network.
nnConcat = neuralNetwork({concatLayer, concatSumLayer});
nnConcat.setInputSize([n0 1]);

% The expected output is 0, because the computation paths cancel each
% other. The cancelation is only observed for point and zonotope
% propagation, not for intervals.
% Check intervals.
Iy = nnConcat.evaluate(Ix);
assert(Iy.inf == -2*n0 & Iy.sup == 2*n0);
% Check points.
ys = nnConcat.evaluate(xs);
assert(all(withinTol(ys,0,1e-7),'all'));
% Check zonotopes.
nnConcat.prepareForZonoBatchEval(xs,options);
[cy,Gy] = nnConcat.evaluateZonotopeBatch( ...
    zeros([n0 bSz]),repmat(eye(n0,q0),1,1,bSz));
assert(all(cy == 0,'all') && all(Gy == 0,'all'));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
