function res = create_nn_unittests()
% create_nn_unittests - creates a reference output for later comparison
%
%
% Syntax:
%    res = create_nn_unittests()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_nn_polyZonotope_adaptive_output_ref

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% SET PARAMETERS

options = struct();
options.nn.bound_approx = true;
options.nn.num_generators = 100;
options.nn.activation = "sigmoid";
options.nn.poly_method = 'taylor';
options.nn.neurons = [2, 5, 5, 2];

name = sprintf( ...
    "./cora/models/Cora/nn/unitTests/model_test_nn_taylm_%s_%s.mat", ...
    options.nn.poly_method, ...
    options.nn.activation ...
    );

% INPUT SET

c = [4; 4];
G = [2, 1, 2; 0, 2, 2];
E = [1, 0, 3; 0, 1, 1];
GI = [];
% input_ref = polyZonotope(c,G,GI,E);
input_ref = taylm(interval([-1;0],[9;8]));

% CREATE NETWORK

neurons = options.nn.neurons;

W = cell(length(neurons)-1, 1);
b = cell(length(neurons)-1, 1);

scale = 1;
for i = 1:length(neurons) - 1
    W{i} = rand(neurons(i+1), neurons(i)) * scale - scale / 2;
    b{i} = rand(neurons(i+1), 1) * scale - scale / 2;
end

layers = cell(length(W)*2, 1);
for i = 1:length(W)
    layers{2*i-1} = nnLinearLayer(W{i}, b{i});

    layer = nnActivationLayer.instantiateFromString(options.nn.activation);
    layers{2*i} = layer;
end
nn_cora = neuralNetwork(layers);

% EVALUATE

output_ref = nn_cora.evaluate(input_ref, options);

% SAVING

save(name, "W", "b", "input_ref", "options", "output_ref")

res = true;

end

% ------------------------------ END OF CODE ------------------------------
