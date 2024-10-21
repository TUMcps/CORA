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
% See also: none

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   14-December-2022 (added variable input tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check basic example
W = rand(4,3); 
b = rand(4,1);
layer = nnLinearLayer(W, b);

assert(compareMatrices(W, layer.W))
assert(compareMatrices(b, layer.b))

% check variable input
layer = nnLinearLayer(W);
assert(sum(layer.b) == 0)

name = "TestLayer";
layer = nnLinearLayer(W, b, name);
assert(layer.name == name)

% wrong input
assertThrowsAs(@nnLinearLayer,'MATLAB:minrhs');

% dimension missmatch
W = rand(4,3); 
b = rand(10,1);
assertThrowsAs(@nnLinearLayer,'CORA:wrongInputInConstructor',W,b);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
