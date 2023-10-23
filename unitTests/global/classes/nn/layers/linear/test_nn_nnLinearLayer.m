function [res] = test_nn_nnLinearLayer()
% test_nn_nnLinearLayer - tests constructor of nnLinearLayer
%
% Syntax:
%    res = test_nn_nnLinearLayer()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   14-December-2022 (added variable input tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check basic example
W = rand(4,3); 
b = rand(4,1);
layer = nnLinearLayer(W, b);

if ~compareMatrices(W, layer.W) && ~compareMatrices(b, layer.b)
    res = false;
end

% check variable input
layer = nnLinearLayer(W);
if sum(layer.b) ~= 0
    res = false;
end

name = "TestLayer";
layer = nnLinearLayer(W, b, name);
if layer.name ~= name
    res = false;
end

% wrong input
try
    layer = nnLinearLayer();
    % should have thrown an error
    res = false;
catch
end

% dimension missmatch
W = rand(4,3); 
b = rand(10,1);
try
    layer = nnLinearLayer(W, b);
    % should have thrown an error
    res = false;
catch
end


end

% ------------------------------ END OF CODE ------------------------------
