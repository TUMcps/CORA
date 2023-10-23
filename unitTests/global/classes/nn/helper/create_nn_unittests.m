function create_nn_unittests()
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

evParams = struct();
evParams.bound_approx = true;
evParams.num_generators = 100;
evParams.activation = "ReLU";
evParams.poly_method = 'regression';
evParams.neurons = [2, 5, 5, 2];

name = sprintf( ...
    "./cora/models/Cora/nn/unitTests/model_test_nn_polyZonotope_%s_%s.mat", ...
    evParams.poly_method, ...
    evParams.activation ...
    );

% INPUT SET

c = [4; 4];
G = [2, 1, 2; 0, 2, 2];
E = [1, 0, 3; 0, 1, 1];
GI = [];
input_ref = polyZonotope(c, G, GI, E);

% CREATE NETWORK

neurons = evParams.neurons;

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

    layer = nnActivationLayer.instantiateFromString(evParams.activation);
    layers{2*i} = layer;
end
nn_cora = neuralNetwork(layers);

% EVALUATE

output_ref = nn_cora.evaluate(input_ref, evParams);

% SAVING

save(name, "W", "b", "input_ref", "evParams", "output_ref")

end

% ------------------------------ END OF CODE ------------------------------
