function example_nn_open_loop_01()
% example_nn_open_loop_01 - runs an open-loop neural network verification 
%    example
%
%
% Syntax:
%    res = example_nn_open_loop_01()
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
% See also: -

% Author:       Tobias Ladner
% Written:      17-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

seed = randi(1000);
seed = 231;
rng(seed)
fprintf("Seed: %d\n", seed);

% INIT

% init random neural network
neurons = [2, 2];

W = cell(length(neurons)-1, 1);
b = cell(length(neurons)-1, 1);

scale = 4;
for i=1:length(neurons)-1
    W{i} = rand(neurons(i+1), neurons(i)) * scale - scale/2;
    b{i} = rand(neurons(i+1), 1) * scale - scale/2;
end

layers = cell(length(W)*2, 1);
for i=1:length(W)
    layers{2*i-1} = nnLinearLayer(W{i}, b{i});
    layers{2*i} = nnSigmoidLayer();
end
nn_cora = neuralNetwork(layers);

% input set
c = [4;4];
G = [2 1 2; 0 2 2];
expMat = [1 0 3;0 1 1];
Grest = [];
S = polyZonotope(c,G,Grest,expMat);

% evaluate random sample points
points = S.randPoint(1000);
res_points = nn_cora.evaluate(points);

% PLOTTING

% plot input set
figure; hold on;
subplot(1, 2, 1); hold on;
plot(S, [1, 2], "DisplayName", "Input Set")
scatter(points(1, 1:end), points(2, 1:end), '.k', "DisplayName", "Samples")
legend('Location', "best")
title("Input")

% plot output set
subplot(1, 2, 2); hold on;
scatter(res_points(1, 1:end), res_points(2, 1:end), '.k', "DisplayName", "Samples")
legend('Location', "best")
title("Output")
drawnow

% EVALUATE NETWORK
types = {'lin', 'quad', 'cub'};
cmap = turbo(length(types));

evParams = struct();
evParams.bound_approx = false;
evParams.num_generators = 1000;

for i=1:length(types)
    evParams.polynomial_approx = types{i};
    
    % evaluate
    fprintf("Evaluating with '%s'\n", types{i});
    res_new = nn_cora.evaluate(S, evParams);

    % plot
    disp("Plotting ...")
    res_new_plot = reduce(res_new, 'girard', 500);
    plot(res_new_plot, [1, 2], '-', "DisplayName", types{i}, "Color", cmap(i, :))
    drawnow
end
end

%------------- END OF CODE --------------