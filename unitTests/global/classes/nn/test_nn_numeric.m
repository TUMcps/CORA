function res = test_nn_numeric()
% test_nn_numeric - tests the reference output using a numeric input
%
% Syntax:
%    res = test_nn_numeric()
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
% Written:       24-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load W, b, input_ref, output_ref (from previous model)
load('model_test_nn_numeric.mat');

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

res = compareMatrices(output_ref, output_new);
end

% ------------------------------ END OF CODE ------------------------------
