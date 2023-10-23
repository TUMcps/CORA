function [res] = test_nn_nnElementwiseAffineLayer()
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
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       14-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% check basic example
scale = 2;
offset = 3;
layer = nnElementwiseAffineLayer(scale, offset);
if layer.scale ~= scale || layer.offset ~= offset
    res = false;
end

scale = [2;2];
offset = [-4;-4];
layer = nnElementwiseAffineLayer(scale, offset);
if ~compareMatrices(layer.scale, scale) || ...
    ~compareMatrices(layer.offset, offset)
    res = false;
end

% check different dimensions
layer = nnElementwiseAffineLayer([2;3], -4);
layer = nnElementwiseAffineLayer(2, [-1; 2]);

% check variable input
layer = nnElementwiseAffineLayer(scale);
if sum(layer.offset) ~= 0
    res = false;
end

name = "TestLayer";
layer = nnElementwiseAffineLayer(scale, offset, name);
if layer.name ~= name
    res = false;
end

% test wrong inputs
try
    layer = nnElementwiseAffineLayer([2, 2], 1);
    % should have thrown an error
    res = false;
end
try
    layer = nnElementwiseAffineLayer(scale, [-4, 4]);
    % should have thrown an error
    res = false;
end
try
    layer = nnElementwiseAffineLayer([-3;2;1], [-4; 4]);
    % should have thrown an error
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
