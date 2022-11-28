function res = test_nn_zonotope_relu()
% test_nn_zonotope_relu - tests nn with relu activation using zonotopes
%
% Syntax:  
%    res = test_nn_zonotope_relu()
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

% Author:       Tobias Ladner
% Written:      24-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load W, b, input_ref, output_ref (from previous model)
load('model_test_nn_zonotope_relu.mat');

% build the new layer based model
layers = {};
for i = 1:length(W)
    % concat layers
    layers = [; ...
        layers; ...
        {nnLinearLayer(W{i}, b{i})}; ...
        {nnReLULayer()}; ...
        ];
end
nn_new = neuralNetwork(layers);

% calculate output
output_new = nn_new.evaluate(input_ref);

res = is_close(output_ref.Z, output_new.Z);
end

%------------- END OF CODE --------------
