function res = test_nn_interval()
% test_nn_interval - tests nn with sigmoid activation using taylor models 
%
% Syntax:
%    res = test_nn_interval()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       16-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

X = interval([-1;-1], [1;1]);
xs = X.randPoint(1000);

res = true;

for act = ["relu", "sigmoid", "tanh"]

    % create random network
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
        layers{2*i} = nnActivationLayer.instantiateFromString(act);
    end
    nn = neuralNetwork(layers);

    ys = nn.evaluate(xs); % evaluate samples
    Y = nn.evaluate(X); % evaluate input set

    % check point containment
    res = res && all(Y.contains(ys));
end

end

% ------------------------------ END OF CODE ------------------------------
