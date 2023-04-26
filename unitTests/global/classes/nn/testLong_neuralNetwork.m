function res = testLong_neuralNetwork()
% testLong_neuralNetwork - unit test function for the
%                                  neuralNetworkOld class
%
% Syntax:  
%    res = testLong_neuralNetwork()
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
% See also: polygon

% Author:       Niklas Kochdumper
% Written:      17-September-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% test if the computed image encloses all results for random initial points
for i = 1:4
    
    % generate random neural network
    nn = neuralNetworkOld.generateRandom('ActivationFun','sigmoid','NrLayers',i);
    
    % generate random initial set
    X0 = 0.01*zonotope.generateRandom('Dimension',nn.nrOfInputs);
    
    % compute image for the network
    X = evaluate(nn,X0);
    
    % evaluate neural network for random initial points
    points = zeros(nn.nrOfOutputs,1000);
    
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

%------------- END OF CODE --------------