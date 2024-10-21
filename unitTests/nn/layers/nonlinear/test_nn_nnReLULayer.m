function res = test_nn_nnReLULayer()
% test_nn_nnReLULayer - tests the ReLU layer
%
% Syntax:
%    res = test_nn_nnReLULayer()
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
layer = nnReLULayer();

% test values
assert(layer.f(0) == 0);
assert(layer.f(inf) == inf);
assert(layer.f(-1000) == 0);

% test alpha
layer = nnReLULayer();
assert(layer.alpha == 0);

% test name
customName = 'MyLayer';
layer = nnReLULayer(customName);
assert(strcmp(layer.name,customName));

% check evaluate
layer = nnReLULayer();

% check point
x = [1;0;-2];
y = layer.evaluate(x);
assert(all([1;0;0] == y));

% check zonotope
X = zonotope(x,0.01 * eye(3));
Y = layer.evaluate(X);

assert(contains(Y,y));

% gather results
res = true;


% ------------------------------ END OF CODE ------------------------------
