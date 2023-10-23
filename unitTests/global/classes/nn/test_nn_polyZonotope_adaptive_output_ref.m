function res = test_nn_polyZonotope_adaptive_output_ref()
% test_nn_polyZonotope_adaptive_output_ref - test if current output matches
%     a previous reference output
%
%
% Syntax:
%    res = test_nn_polyZonotope_adaptive_output_ref
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: create_nn_unittests

% Authors:       Tobias Ladner
% Written:       24-June-2022
% Last update:   27-April-2023 (added tol to isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];
for activation = ["sigmoid", "tanh", "ReLU"]

    % LOAD SNAPSHOT OUTPUT

    name = sprintf( ...
        "model_test_nn_polyZonotope_%s_%s.mat", ...
        "adaptive", ...
        activation ...
        );
    load(name)

    % CREATE NETWORK

    layers = cell(length(W)*2, 1);
    for i = 1:length(W)
        layers{2*i-1} = nnLinearLayer(W{i}, b{i});

        layer = nnActivationLayer.instantiateFromString(evParams.activation);
        layers{2*i} = layer;
    end
    nn_cora = neuralNetwork(layers);

    % EVALUATE

    output_new = nn_cora.evaluate(input_ref, evParams);
    output_ref = output_ref.replaceId(output_ref.id, (1:length(output_ref.id))');

    % TEST OUTPUT
    resvec(end+1) = isequal(output_ref, output_new,1e-15);

    % figure; hold on;
    % plot(output_ref, [1, 2], 'r')
    % plot(output_new, [1, 2], 'g')
    % title(activation)
    % samples = input_ref.randPoint(5000);
    % pred = nn_cora.evaluate(samples);
    % scatter(pred(1, :), pred(2, :), '.k')
end

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
