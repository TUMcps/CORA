function res = test_nn_conZonotope_ReLU_approx()
% test_nn_conZonotope_ReLU_approx - tests nn with relu activation using
% conZonotopes (approx)
%
% Syntax:
%    res = test_nn_conZonotope_ReLU_approx()
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
model = "model_test_nn_conZonotope_ReLU_approx.mat";
load(model)

% build the new layer based model
layers = {};
for i = 1:length(W)
    % concat layers
    activation_layer = nnActivationLayer.instantiateFromString(options.nn.activation);
    layers = [; ...
        layers; ...
        {nnLinearLayer(W{i}, b{i})}; ...
        {activation_layer}; ...
        ];
end
nn_new = neuralNetwork(layers);

% calculate output
output_new = nn_new.evaluate(input_ref, options);

assert(isequal(output_ref, output_new));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
