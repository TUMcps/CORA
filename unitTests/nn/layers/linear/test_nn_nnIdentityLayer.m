function [res] = test_nn_nnIdentityLayer()
% test_nn_nnIdentityLayer - tests constructor of 
%    nnIdentityLayer
%
% Syntax:
%    res = test_nn_nnIdentityLayer()
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

% simple example
layer = nnIdentityLayer();

% check evaluate

% check point
x = [1;2;3;4];
y = layer.evaluate(x);

assert(all(x == y));

% check zonotope
X = zonotope(x,0.01 * eye(4));
Y = layer.evaluate(X);

assert(isequal(X,Y));

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
