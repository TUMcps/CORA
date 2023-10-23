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

resvec = [];

% simple example
layer = nnIdentityLayer();
resvec(end+1) = true;

% check evaluate

% check point
x = [1;2;3;4];
y = layer.evaluate(x);

resvec(end+1) = all(x == y);

% check zonotope
X = zonotope(x,0.01 * eye(4));
Y = layer.evaluate(X);

resvec(end+1) = isequal(X,Y);

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
