function example_neuralNetwork_evaluate_adaptive()
% example_neuralNetwork_evaluate_adaptive - example for 
%    neural network verification using an unsafe set as specification
%
% Syntax:
%    res = example_neuralNetwork_evaluate_adaptive()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Authors:       Tobias Ladner
% Written:       02-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

seed = 17;
rng(seed)
fprintf("Seed: %d\n", seed);

% Create neural network ---------------------------------------------------

neurons = [2, 5, 5, 2];

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

% Create input set --------------------------------------------------------

c = [4;4];
G = [2 1 2; 0 2 2];
E = [1 0 3;0 1 1];
GI = [];
X = polyZonotope(c,G,GI,E);

% propagate samples
xs = X.randPoint(5000);
ys = nn_cora.evaluate(xs);

% Plot --------------------------------------------------------------------

% only plot if called directly
doPlot = length(dbstack) == 1;

if doPlot
    % plot input
    figure;
    subplot(1, 2, 1); hold on;
    plot(X, [1, 2], "DisplayName", "Input Set")
    scatter(xs(1, 1:end), xs(2, 1:end), '.k', "DisplayName", "Samples")
    legend('Location', "best")
    title("Input")
    
    % plot output
    subplot(1, 2, 2); hold on;
    legend('Location', "best")
    title("Output")
    scatter(ys(1, 1:end), ys(2, 1:end), '.k', "DisplayName", "Samples")
    drawnow
end

% Run Adaptive Evaluation -------------------------------------------------

disp("Adaptive")
steps = 7;
cmap = turbo(steps);

evParams = struct();
evParams.bound_approx = false;
evParams.reuse_bounds = true; % change to 'false' for more accurate result
evParams.num_generators = 10000;
evParams.poly_method = "regression";
evParams.max_bounds = 5;

for i=1:steps
    fprintf("Step %i\n", i);
    Y = nn_cora.evaluate(X, evParams);
    
    % plot
    if doPlot
    % reduce for plotting
        Y = reduce(Y, 'girard', 500);
        
        % plot order used per neuron in legend
        refs = [];
        refinable_layers = nn_cora.getRefinableLayers();
        for j=1:length(refinable_layers)
            refs = [refs, ' [', num2str((refinable_layers{j}.order)'), ']'];
        end
    
        disp("Plotting ...")
        plot(Y, [1, 2], '-', "DisplayName", sprintf(['Order Pattern:', refs]), "Color", cmap(i, :))
        drawnow
    end
    
    % adaptively refine network
    if i < steps
        disp("Refining ...")
        nn_cora.refine(evParams.max_bounds, "layer", "both", X.randPoint(1), true, [])
    end
end

end

% ------------------------------ END OF CODE ------------------------------
