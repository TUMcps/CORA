function res = test_nn_neuralNetwork_evaluate()
% test_nn_neuralNetwork_evaluate - unit test function for 
%     neuralNetwork/evaluate: sanity check of points
%
% Syntax:
%    res = test_nn_neuralNetwork_evaluate()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   27-April-2023 (TL, removed neuralNetworkOld)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% Test 1: Computed image encloses all results for random initial points
for i = 1:4
    
    % generate random neural network
    nn = neuralNetwork.generateRandom('ActivationFun','sigmoid','NrLayers',i);
    
    % generate random initial set
    X0 = 0.01*zonotope.generateRandom('Dimension',nn.neurons_in);
    
    % compute image for the network
    Y = evaluate(nn,X0);
    
    % evaluate neural network for random initial points
    xs = randPoint(X0,100);
    ys = nn.evaluate(xs);
    
    % check if all points are inside the computed image
    resvec(end+1) = all(contains(Y,ys));
end

% Test 2: Propagate different sets

nn = neuralNetwork.generateRandom('ActivationFun','sigmoid','NrInputs',2);

I = interval.generateRandom("Dimension",2);
nn.evaluate(I);
resvec(end+1) = true;

Z = zonotope(I);
nn.evaluate(Z);
resvec(end+1) = true;

pZ = polyZonotope(I);
nn.evaluate(pZ);
resvec(end+1) = true;

pZ = polyZonotope(Z.c,[],Z.G);
nn.evaluate(pZ);
resvec(end+1) = true;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
