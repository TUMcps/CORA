function res = test_nn_nnActLayerFromHandle()
% test_nn_nnActLayerFromHandle - tests the nnActLayerFromHandle
%
% Syntax:
%    res = test_nn_nnActLayerFromHandle()
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
% Written:       06-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

funs = {@sin, @cos};

for i=1:numel(funs)
    % init layer
    layer = nnActLayerFromHandle(funs{i});
    
    % check evaluate
    
    % check point
    x = [0;0];
    y = layer.evaluate(x);
    assertLoop(all(layer.f(x) == y),i);
    
    % check zonotope
    X = zonotope(x,pi/3 * eye(2));
    Y = layer.evaluate(X);
    % check containment of points
    xs = X.randPoint(100);
    ys = layer.evaluate(xs);
    assertLoop(contains(Y,ys),i);
end

% gather results
res = true;


% ------------------------------ END OF CODE ------------------------------
