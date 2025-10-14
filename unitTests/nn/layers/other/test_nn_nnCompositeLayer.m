function [res] = test_nn_nnCompositeLayer()
% test_nn_nnCompositeLayer - tests constructor of nnCompositeLayer.
%
% Syntax:
%    res = test_nn_nnCompositeLayer()
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

% Authors:       Lukas Koller
% Written:       27-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Instantiate layers.
W1 = rand(4,3); 
b1 = rand(4,1);
layer1 = nnLinearLayer(W1, b1);

W2 = rand(5,3); 
b2 = rand(5,1);
layer2 = nnLinearLayer(W2, b2);

W3 = rand(4,5); 
b3 = rand(4,1);
layer3 = nnLinearLayer(W3, b3);

layers = {
    {layer1, nnReLULayer};
    {layer2, nnTanhLayer, layer3}
};

% Check layers.
compLayer = nnCompositeLayer(layers,'add');
assert(isequal(compLayer.layers,layers));

% Check aggregation.
compLayer = nnCompositeLayer(layers,'add');
assert(strcmp(compLayer.aggregation,'add'));

compLayer = nnCompositeLayer(layers,'concat');
assert(strcmp(compLayer.aggregation,'concat'));

% Wrong aggregation.
assertThrowsAs(@nnCompositeLayer, ...
    'CORA:wrongInputInConstructor',layers,'minus');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
