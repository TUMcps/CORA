function res = test_nn_nnCompositeLayer_concat_evaluate()
% test_nn_nnCompositeLayer_concat_evaluate - evaluation of a 'concat'
%    nnCompositeLayer with ReLU paths: output equals [f(x); g(x)] for point
%    evaluation and soundly encloses it for interval and zonotope-batch
%    propagation
%
% Syntax:
%    res = test_nn_nnCompositeLayer_concat_evaluate()
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
% See also: nnCompositeLayer, nnHelper.buildJointNetwork

% Authors:       Benedikt Kellner
% Written:       03-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng(7);

nin = 4;
mf = 3;   % output dim of path f
mg = 2;   % output dim of path g

% Two independent ReLU sub-networks sharing the same input.
fLayers = {nnLinearLayer(rand(6,nin)-0.5, rand(6,1)-0.5), nnReLULayer(), ...
           nnLinearLayer(rand(mf,6)-0.5, rand(mf,1)-0.5)};
gLayers = {nnLinearLayer(rand(5,nin)-0.5, rand(5,1)-0.5), nnReLULayer(), ...
           nnLinearLayer(rand(mg,5)-0.5, rand(mg,1)-0.5)};

% Reference networks (same input x).
nn_f = neuralNetwork(fLayers); nn_f.setInputSize([nin 1]);
nn_g = neuralNetwork(gLayers); nn_g.setInputSize([nin 1]);

% Concatenated composite network.
compLayer = nnCompositeLayer({fLayers; gLayers}, 'concat');
nn = neuralNetwork({compLayer});
nn.setInputSize([nin 1]);

% Output size bookkeeping.
assert(isequal(compLayer.getOutputSize([nin 1]), [mf+mg 1]), ...
    'concat output size should be mf+mg');

% --- 1. Point evaluation: exact equivalence to [f(x); g(x)] ---
X = 2*rand(nin, 100) - 1;
Y = nn.evaluate(X);
Yref = [nn_f.evaluate(X); nn_g.evaluate(X)];
assert(isequal(size(Y), [mf+mg, 100]), 'concat point output has wrong size');
assert(all(abs(Y - Yref) < 1e-9, 'all'), ...
    'concat point output != [f(x); g(x)]');

% --- 2. Interval propagation: sound enclosure of point outputs ---
Ix = interval(-ones(nin,1), ones(nin,1));
Iy = nn.evaluate(Ix);
assert(dim(Iy) == mf+mg, 'concat interval output has wrong dim');
% f- and g-blocks must each enclose their sub-network's interval image.
Iyf = nn_f.evaluate(Ix);
Iyg = nn_g.evaluate(Ix);
Iref = [Iyf; Iyg];
assert(all(infimum(Iy) <= infimum(Iref) + 1e-7), 'concat interval lower not sound');
assert(all(supremum(Iy) >= supremum(Iref) - 1e-7), 'concat interval upper not sound');
% And it must enclose sampled point outputs.
Ys = nn.evaluate(2*rand(nin, 200) - 1);
assert(all(Ys >= infimum(Iy) - 1e-7, 'all') && ...
       all(Ys <= supremum(Iy) + 1e-7, 'all'), ...
       'concat interval does not enclose sampled outputs');

% --- 3. Zonotope-batch propagation: sound enclosure of point outputs ---
q0 = nin;
options.nn.train.num_init_gens = nin;
% Representative points to set up batch evaluation state.
xs = 2*rand(nin, 100) - 1;
nn.prepareForZonoBatchEval(xs, options);
bSz = 3;
c = zeros(nin, bSz);                       % zonotope centred at origin
G = repmat(eye(nin, q0), 1, 1, bSz);       % identity generators -> box [-1,1]^nin
[cy, Gy] = nn.evaluateZonotopeBatch(c, G);

assert(size(cy,1) == mf+mg, 'concat zonotope center has wrong dim');
assert(size(Gy,1) == mf+mg, 'concat zonotope generators have wrong dim');

% Interval over-approximation of the output zonotope (first batch element).
rad = sum(abs(Gy(:,:,1)), 2);
yl = cy(:,1) - rad;
yu = cy(:,1) + rad;
% Sample points from the input box and check enclosure of point outputs.
Xz = 2*rand(nin, 500) - 1;
Yz = nn.evaluate(Xz);
assert(all(Yz >= yl - 1e-6, 'all') && all(Yz <= yu + 1e-6, 'all'), ...
    'concat zonotope does not enclose sampled outputs');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
