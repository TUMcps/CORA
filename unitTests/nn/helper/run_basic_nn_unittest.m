function res = run_basic_nn_unittest(optionsCustom)
% run_basic_nn_unittest - test basic functionality of nn evaluation
%
%
% Syntax:
%    res = run_basic_nn_unittest(optionsCustom)
%
% Inputs:
%    optionsCustom - struct (see neuralNetwork/evaluate)
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate, neuralNetwork/refine

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   29-November-2022 (numRefine)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

seed = 17;
rng(seed);

% SET PARAMETERS

options = struct();
options.nn.bound_approx = true;
options.nn.num_generators = 1000;
options.nn.reuse_bounds = false;

options.nn.neurons = [2, 5, 5, 2];
options.nn.activation = "sigmoid";

options.nn.refine = false;
options.nn.numRefine = 1;

options.nn.max_bounds = 5;
options.nn.refine_type = "layer";
options.nn.refine_method = "approx_error";

for field = fieldnames(optionsCustom.nn)
    options.nn.(field{1}) = optionsCustom.nn.(field{1});
end

% INPUT SET

c = [4; 4];
G = [2, 1, 2; 0, 2, 2];
E = [1, 0, 3; 0, 1, 1];
GI = [];
S_in = polyZonotope(c, G, GI, E);

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
    if isfield(options.nn, 'order')
        layer.order = options.nn.order;
    end
    layers{2*i} = layer;
end
nn_cora = neuralNetwork(layers);

% RUN EVALUATE

for i=1:options.nn.numRefine
    S_out = nn_cora.evaluate(S_in, options);
    
    if options.nn.refine
        nn_cora.refine(options.nn.max_bounds, options.nn.refine_type, options.nn.refine_method, S_in.randPoint(1))
        S_out = nn_cora.evaluate(S_in, options);
    end
end

% TEST FOR POINTS

P_in = [S_in.randPoint(100), S_in.randPoint(50, 'extreme')];
P_out = nn_cora.evaluate(P_in, options);

assert(all(zonotope(S_out).contains(P_out)));

res = true;

% figure; hold on;
% plot(S_out);
% scatter(P_out(1, :), P_out(2, :), '.k');
% drawnow;

end

% ------------------------------ END OF CODE ------------------------------
