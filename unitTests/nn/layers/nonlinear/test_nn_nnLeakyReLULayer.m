function res = test_nn_nnLeakyReLULayer()
% test_nn_nnLeakyReLULayer - tests the leaky ReLU layer
%
% Syntax:
%    res = test_nn_nnLeakyReLULayer()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       02-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init layer
layer = nnLeakyReLULayer();

% test values
assert(layer.f(0) == 0);
assert(layer.f(inf) == inf);
assert(layer.f(-inf) == -inf);

% test alpha
alpha = -0.1;
layer = nnLeakyReLULayer(alpha);
assert(layer.alpha == alpha);

% test name
customName = 'MyLayer';
layer = nnLeakyReLULayer(0.01, customName);
assert(strcmp(layer.name,customName));

% check evaluate
alpha = 0.01;
layer = nnLeakyReLULayer(alpha);

% check point
x = [1;0;-2];
y = layer.evaluate(x);
assert(all([1;0;-2*alpha] == y));

% check zonotope
X = zonotope(x,0.01 * eye(3));
Y = layer.evaluate(X);

assert(contains(Y,y));

% gather results
res = true;


% ------------------------------ END OF CODE ------------------------------
