function res = test_nn_nnRootLayer()
% test_nn_nnRootLayer - tests the square root layer
%
% Syntax:
%    res = test_nn_nnRootLayer()
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
% Written:       02-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init layer
layer = nnRootLayer();

% test name
customName = 'MyLayer';
layer = nnRootLayer(customName);
assert(strcmp(layer.name,customName));

% check evaluate

% check point
x = [1;2;3;4];
y = layer.evaluate(x);
assert(all(layer.f(x) == y));

% check zonotope
X = zonotope(x,0.01 * eye(4));
Y = layer.evaluate(X);

assert(contains(Y,y));

% gather results
res = true;


% ------------------------------ END OF CODE ------------------------------
