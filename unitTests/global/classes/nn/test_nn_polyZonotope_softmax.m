function res = test_nn_polyZonotope_softmax()
% test_nn_polyZonotope_softmax - test softmax layer
%
%
% Syntax:
%    res = test_nn_polyZonotope_softmax
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
% Written:       17-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.num_generators = 100;

% INPUT SET

S_in = polyZonotope.generateRandom('Dimension',10,'NrGenerators',15);

% CREATE NETWORK

layers = {nnSoftmaxLayer()};
nn_cora = neuralNetwork(layers);

% RUN EVALUATE

S_out = nn_cora.evaluate(S_in, evParams);

% test sensitivity
nn_cora.calcSensitivity(S_in.c);

% TEST FOR POINTS

P_in = [S_in.randPoint(100), S_in.randPoint(50, 'extreme')];
P_out = nn_cora.evaluate(P_in, evParams);

resvec = contains(zonotope(S_out),P_out);
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
