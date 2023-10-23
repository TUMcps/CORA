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

resvec = [];

% init layer
layer = nnReLULayer();

% test values
resvec(end+1) = layer.f(0) == 0;
resvec(end+1) = layer.f(inf) == inf;
resvec(end+1) = layer.f(-1000) == 0;

% test alpha
layer = nnReLULayer();
resvec(end+1) = layer.alpha == 0;

% test name
customName = 'MyLayer';
layer = nnReLULayer(customName);
resvec(end+1) = strcmp(layer.name,customName);

% check evaluate
layer = nnReLULayer();

% check point
x = [1;0;-2];
y = layer.evaluate(x);
resvec(end+1) = all([1;0;0] == y);

% check zonotope
X = zonotope(x,0.01 * eye(3));
Y = layer.evaluate(X);

resvec(end+1) = contains(Y,y);

% gather results
res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
