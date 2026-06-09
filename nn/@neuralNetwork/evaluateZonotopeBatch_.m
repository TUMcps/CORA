function [c,G,a] = evaluateZonotopeBatch_(nn,c,G,options,idxLayer)
% evaluateZonotopeBatch_ - evaluate neural network for a batch of zonotopes
%   without setting default options.
%
% Syntax:
%    [c, G] = nn.evaluateZonotopeBatch_(c, G, options, idxLayer)
%
% Inputs:
%    c, G - batch of zonotope; [n,q+1,b] = size([c G]),
%       where n is the number of dims, q the number of generators, and b the batch size
%    options - parameter for neural network evaluation
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    c, G - batch of output sets
%    a - struct, aggreagation result from options.nn.neuron_aggregation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate, neuralNetwork/prepareForZonoBatchEval

% Authors:       Lukas Koller
% Written:       02-August-2023
% Last update:   08-August-2023 (moved code to layers)
%                22-February-2024 (merged options.nn, moved input storage handling from layer to network)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Initialize the neuron aggregation result.
a = struct();
% Check if there is a neuron aggregation function.
if isfield(options.nn,'neuron_aggregation_fun') && nargout == 3
    doNeuronAggregation = true; 
    % Extract the neuron aggregation function.
    neurAggFun = options.nn.neuron_aggregation_fun;
    % Re-set the neuron aggregation results in the options; we pass it 
    % using the options to nnCompositeLayer/evaluteZonotopeBatch.
    options.nn.neuron_aggregation_result = a;
else
    % There is no neuron aggregation.
    doNeuronAggregation = false; 
end

for i=idxLayer
    % Obtain the i-th layer.
    layeri = nn.layers{i};

    % Store input for neuron-splitting or backpropgation.
    if options.nn.train.backprop || ...
            (isa(layeri,'nnActivationLayer') && options.nn.backprop_without_weight_update)
        layeri.backprop.store.inc = c;
        layeri.backprop.store.inG = G;
    end

    % Save pre-activation input for aggregation.
    if doNeuronAggregation
        cin = c; 
        Gin = G;
    end

    if isa(layeri,'nnCompositeLayer') && doNeuronAggregation
        % Update the neuron aggregation results in the options for
        % nnCompositeLayer/evaluteZonotopeBatch.
        options.nn.neuron_aggregation_result = a;
        % Compute the result of the i-th layer.
        [c,G,a] = layeri.evaluateZonotopeBatch(c,G,options);
    else
        % Compute the result of the i-th layer.
        [c,G] = layeri.evaluateZonotopeBatch(c,G,options);
    end

    if doNeuronAggregation
        % Call the neuron aggregation function (after layer evaluation,
        % passing both pre-activation input and post-activation output).
        a = neurAggFun(a,layeri,cin,Gin,c,G);
    end
end

if doNeuronAggregation
    % Clear the neuron aggregation result.
    options = rmfield(options.nn,'neuron_aggregation_result');
end

end

% ------------------------------ END OF CODE ------------------------------
