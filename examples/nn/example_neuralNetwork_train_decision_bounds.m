function res = example_neuralNetwork_train_decision_bounds()
% example_neuralNetwork_train_decision_bounds - example for training a neural network.
%
% Syntax:
%    res = example_neuralNetwork_train_decision_bounds()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%    [1] Koller, L. et al. Set-based Training for Neural Network
%           Verification. TMLR 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       12-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Reset the random number generator.
rng(42)

% 1. Generate Data. -------------------------------------------------------

% Define the input space.
n0 = 2; % Specify the number of input dimensions
% Specify the bounds of the input space.
l = [0; 0];
u = [1; 1];

% Specify perturbation radius.
epsilon = 0.05;

% Specify the number of data points.
N = 20;

% Sample non-intersecting data points.
combs = combvec(l(1)+2*epsilon:4*epsilon:u(1)-2*epsilon,...
    l(2)+2*epsilon:4*epsilon:u(2)-2*epsilon);
ids = randperm(length(combs),N);
xs = combs(:,ids) + 2*epsilon*(rand(n0,N) - 1/2);

% Define ground truth.
ts = zeros(n0,N);
ts(sub2ind(size(ts),randi([1 2],1,N),1:N)) = 1;
% Extract target labels.
[~,ls] = max(ts);

% 2. Training (point-based & set-based). ----------------------------------

% Construct a neural network.
layers = aux_initNeuralNetworkLayers();
nn1 = neuralNetwork(layers);
% Create a copy of the initialized neural network.
nn2 = nn1.copyNeuralNetwork();

% Specify training options.
options.nn.train = struct( ...
    'optim',nnAdamOptimizer(1e-3),... 
    'max_epoch',500,...
    'mini_batch_size',10,...
    'loss','softmax+log',...
    'shuffle_data','every_epoch',...
    'noise',epsilon,...
    'warm_up',100,...
    'ramp_up',450, ...
    'tau',0.01,...
    'init_gens','l_inf',...
    'num_init_gens',inf,...
    'num_approx_err',inf,...
    'volume_heuritsic','f-radius',...
    'zonotope_weight_update','sum',...
    'exact_backprop',true ...
);

% Apply point-based training.
options1 = options;
options1.nn.train.method = 'point';
loss1 = nn1.train(xs,ts,zeros([2 0]),[],options1);

% % Plot the loss.
% figure; hold on;
% title(sprintf("Loss ('%s')",options1.nn.train.method))
% plot(1:length(loss1.center),smooth(loss1.center), ...
%     'DisplayName','Training Loss');
% legend

% Compute and print the accuracy.
fprintf('Accuracy (%s): %.2f\n',options1.nn.train.method, ...
    aux_computeAcc(nn1,xs,ls));

% Apply set-based training.
options2 = options;
options2.nn.train.method = 'set';
loss2 = nn2.train(xs,ts,zeros([n0 0]),[],options2);

% % Plot the loss.
% figure; hold on;
% title(sprintf("Loss ('%s')",options2.nn.train.method))
% plot(1:length(loss2.center),smooth(loss2.center), ...
%     'DisplayName','Training Loss (Accuracy)');
% yyaxis right
% plot(1:length(loss2.vol),smooth(loss2.vol), ...
%     'DisplayName','Training Loss (Robustness)');
% legend

% Compute and print the accuracy.
fprintf('Accuracy (%s): %.2f\n',options2.nn.train.method, ...
    aux_computeAcc(nn2,xs,ls));

% Visualize the Decision Bounds. ------------------------------------------

% We approximate the decision bounds with evenly spaced samples. Specify 
% the number of samples for the mesh grid.
numSamples = 100;
% Compute contour matrices for networks.
[xs_1,xs_2,ys1_] = aux_computeNNDecisionBounds(nn1,l,u,numSamples);
[~,~,ys2_] = aux_computeNNDecisionBounds(nn2,l,u,numSamples);

% Plot the decision bounds.
figure;
subplot(1,2,1); hold on;
title(sprintf("Training Method '%s'",options2.nn.train.method))
aux_plotDecisionBounds(xs,ls,l,u,epsilon,xs_1,xs_2,ys1_,[]);

subplot(1,2,2); hold on;
title(sprintf("Training Method '%s'",options2.nn.train.method))
aux_plotDecisionBounds(xs,ls,l,u,epsilon,xs_1,xs_2,ys2_,ys1_);

res = true;

end


% Auxiliary functions -----------------------------------------------------

function layers = aux_initNeuralNetworkLayers()
    % Construct the layers of the neural network.
    dltLayers = [
        featureInputLayer(2)
        fullyConnectedLayer(200)
        reluLayer
        fullyConnectedLayer(200)
        reluLayer
        fullyConnectedLayer(200)
        reluLayer
        fullyConnectedLayer(2)
    ];
    % Initialize the weights and biases using the Deep Learning Toolbox.
    nn_dlt = dlnetwork(dltLayers);
    % Convert the initialized neural network to CORA.
    nn = neuralNetwork.convertDLToolboxNetwork(num2cell(nn_dlt.Layers));
    % Extract the initialized layers.
    layers = nn.layers;

end

function acc = aux_computeAcc(nn,xs,ls)
    % Compute the accuracy [in percent] of a neural network w.r.t. given 
    % data points xs and target labels ls.

    % Obtain the number of data points.
    N = size(xs,2);
    % Compute the outputs of the neural network.
    y = nn.evaluate(xs);
    % Obtain the class predictions.
    [~,ks] = max(y);
    % Compute the accuracy.
    acc = sum(ls == ks)/N*100;
end

function [xs_1,xs_2,ys_] = aux_computeNNDecisionBounds(nn,l,u,N)
    % Approximate the decision bounds of a given neural network using 
    % evenly spaced samples.

    % Specify a batch size.
    bs = 128;

    % Sample the input space evenly.
    xs_1 = linspace(l(1),u(1),N);
    xs_2 = linspace(l(2),u(2),N);

    % Compute mesh grid.
    xs = combvec(xs_1,xs_2);
    ys = zeros(2,size(xs,2));

    % Compute output for mesh.
    for i=1:bs:size(xs,2)
        % Compute the batch indices.
        batchIdx = i:min(i+bs,size(xs,2));
        % Compute the outputs for the current batch.
        ys(:,batchIdx) = nn.evaluate(xs(:,batchIdx));
    end
    % Convert outputs to contour matrix.
    ys_ = reshape([1 -1]*(softmax(ys)),[N N])';
end

function ys_ = aux_plotDecisionBounds(xs,ls,l,u,epsilon, ...
    xs_1,xs_2,ys_,ys_dashed)
    % Plot a given decision boundary.

    % Specify colors.
    color1 = CORAcolor('CORA:highlight1');
    color2 = CORAcolor('CORA:reachSet');

    % Specify the color map.
    nc = 20; % Specify number of colors.
    cs1 = [linspace(color1(1),1,nc); ...
        linspace(color1(2),1,nc); ...
        linspace(color1(3),1,nc)]';
    cs2 = [linspace(color2(1),1,nc); ...
        linspace(color2(2),1,nc); ...
        linspace(color2(3),1,nc)]';

    % Plot contour.
    contourf(xs_1,xs_2,ys_,100,'EdgeColor','none','FaceAlpha',0.5);
    colormap([cs2(1:end-1,:); flip(cs1)])
    colorbar

    % Plot dashed contour.
    if ~isempty(ys_dashed)
        contour(xs_1,xs_2,ys_dashed,[0 0],'--k');
    end
    % Plot data points.
    for i=1:size(xs,2)
        if ls(i) == 1
            scatter(xs(1,i),xs(2,i),50,'filled','MarkerFaceColor',color1)
        else
            scatter(xs(1,i),xs(2,i),50,'filled','MarkerFaceColor',color2)
        end
        plot(xs(:,i) + epsilon*interval([-1;-1],[1;1]),[1 2],'-k')
    end
    xlim([l(1) u(1)])
    ylim([l(2) u(2)])
end

% ------------------------------ END OF CODE ------------------------------
