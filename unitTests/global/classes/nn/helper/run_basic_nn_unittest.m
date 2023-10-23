function res = run_basic_nn_unittest(evParamsCustom)
% run_basic_nn_unittest - test basic functionality of nn evaluation
%
%
% Syntax:
%    res = run_basic_nn_unittest(evParamsCustom)
%
% Inputs:
%    evParamsCustom - struct (see neuralNetwork/evaluate)
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

evParams = struct();
evParams.bound_approx = true;
evParams.num_generators = 1000;
evParams.reuse_bounds = false;

evParams.neurons = [2, 5, 5, 2];
evParams.activation = "sigmoid";

evParams.refine = false;
evParams.numRefine = 1;

evParams.max_bounds = 5;
evParams.refine_type = "layer";
evParams.refine_method = "approx_error";

for field = fieldnames(evParamsCustom)
    evParams.(field{1}) = evParamsCustom.(field{1});
end

% INPUT SET

c = [4; 4];
G = [2, 1, 2; 0, 2, 2];
E = [1, 0, 3; 0, 1, 1];
GI = [];
S_in = polyZonotope(c, G, GI, E);

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
    if isfield(evParams, 'order')
        layer.order = evParams.order;
    end
    layers{2*i} = layer;
end
nn_cora = neuralNetwork(layers);

% RUN EVALUATE

for i=1:evParams.numRefine
    S_out = nn_cora.evaluate(S_in, evParams);
    
    if evParams.refine
        nn_cora.refine(evParams.max_bounds, evParams.refine_type, evParams.refine_method, S_in.randPoint(1))
        S_out = nn_cora.evaluate(S_in, evParams);
    end
end

% TEST FOR POINTS

P_in = [S_in.randPoint(100), S_in.randPoint(50, 'extreme')];
P_out = nn_cora.evaluate(P_in, evParams);

res = true;
res = res && all(zonotope(S_out).contains(P_out));

% figure; hold on;
% plot(S_out);
% scatter(P_out(1, :), P_out(2, :), '.k');
% drawnow;

end

% ------------------------------ END OF CODE ------------------------------
