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

% test if the computed image encloses all results for random initial points
for i = 1:4
    
    % generate random neural network
    nn = neuralNetwork.generateRandom('ActivationFun','sigmoid','NrLayers',i);
    
    % generate random initial set
    X0 = 0.01*zonotope.generateRandom('Dimension',nn.neurons_in);
    
    % compute image for the network
    X = evaluate(nn,X0);
    
    % evaluate neural network for random initial points
    points = zeros(nn.neurons_out,100);
    
    for j = 1:size(points,2)
       p = randPoint(X0);
       points(:,j) = evaluate(nn,p);
    end
    
    % check if all points are inside the computed image
    if ~contains(X,points)
        throw(CORAerror('CORA:testFailed'));
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
