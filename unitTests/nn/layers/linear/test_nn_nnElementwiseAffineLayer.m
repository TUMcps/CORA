function res = test_nn_nnElementwiseAffineLayer()
% test_nn_nnElementwiseAffineLayer - tests constructor of 
%    nnElementwiseAffineLayer
%
% Syntax:
%    res = test_nn_nnElementwiseAffineLayer()
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
% Written:       14-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check basic example
scale = 2;
offset = 3;
layer = nnElementwiseAffineLayer(scale, offset);
assert(layer.scale == scale)
assert(layer.offset == offset)

scale = [2;2];
offset = [-4;-4];
layer = nnElementwiseAffineLayer(scale, offset);
assert(compareMatrices(layer.scale, scale))
assert(compareMatrices(layer.offset, offset))

% check different dimensions
layer = nnElementwiseAffineLayer([2;3], -4);
layer = nnElementwiseAffineLayer(2, [-1; 2]);

% check variable input
layer = nnElementwiseAffineLayer(scale);
assert(sum(layer.offset) == 0)

name = "TestLayer";
layer = nnElementwiseAffineLayer(scale, offset, name);
assert(layer.name == name)

% test wrong inputs
assertThrowsAs(@nnElementwiseAffineLayer,'CORA:wrongInputInConstructor',[2,2],1);
assertThrowsAs(@nnElementwiseAffineLayer,'CORA:wrongInputInConstructor',scale,[-4,4]);
assertThrowsAs(@nnElementwiseAffineLayer,'CORA:wrongInputInConstructor',[-3;2;1],[-4;4]);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
