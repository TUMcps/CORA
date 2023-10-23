function [res] = test_nn_nnAvgPool2DLayer()
% test_nn_nnAvgPool2DLayer - tests constructor of 
%    nnAvgPool2DLayer
%
% Syntax:
%    res = test_nn_nnAvgPool2DLayer()
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

% simple example
layer = nnAvgPool2DLayer([2,2]);
resvec(end+1) = isequal(layer.W,0.25 * ones(2,2));
resvec(end+1) = isequal(layer.stride,[2 2]);

layer = nnAvgPool2DLayer([1,1]);
resvec(end+1) = isequal(layer.W,1);
resvec(end+1) = isequal(layer.stride,[1 1]);

% check evaluate
layer = nnAvgPool2DLayer([2,2]);
nn = neuralNetwork({layer});
n = 4;
nn.setInputSize([n,n,1]);

% check point
x = reshape(eye(n),[],1);
y = nn.evaluate(x);
y_true = [0.5 0 0 0.5]';

resvec(end+1) = all(y == y_true);

% check zonotope
X = zonotope(x,0.01 * eye(n*n));
Y = nn.evaluate(X);

resvec(end+1) = contains(Y,y);

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
