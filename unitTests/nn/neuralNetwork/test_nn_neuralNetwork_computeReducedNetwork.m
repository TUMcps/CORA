function [res] = test_nn_neuralNetwork_computeReducedNetwork()
% test_nn_neuralNetwork_computeReducedNetwork - tests the 
%    computeReducedNetwork function 
%    
%
% Syntax:
%    res = test_nn_neuralNetwork_computeReducedNetwork()
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
% Written:       09-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

seed = 1;
rng(seed);

% create simple network
W1 = rand(10, 2) *2 -1;
b1 = rand(10, 1);

W2 = rand(3, 10) *2 -1;
b2 = rand(3, 1); 

% change bias to have near 0/1 neurons
b1(1:3) = 1000;
b1(end-3:end) = -1000;

nn = neuralNetwork({ ...
    nnLinearLayer(W1, b1);
    nnSigmoidLayer();
    nnLinearLayer(W2, b2);
    nnSigmoidLayer();
});

% get input set
X = zonotope.generateRandom("Dimension", 2);
X = polyZonotope(X);

% sample points
xs = [X.randPoint(100), X.randPoint(100, 'extreme')];
ys = nn.evaluate(xs);

% reduce network
[nn_red, Y_red] = nn.computeReducedNetwork(X);

% check if points are contained
assert(all(zonotope(Y_red).contains(ys)));

% double check with verbose output, should not throw an error
[nn_red, Y_red] = nn.computeReducedNetwork(X, "Verbose",true);

% check if points are contained
assert(all(zonotope(Y_red).contains(ys)));

% reduction rate
[nn_red, Y_red] = nn.computeReducedNetwork(X, "ReductionRate",0.2);

% check if points are contained
assert(all(zonotope(Y_red).contains(ys)));

% input compression
X = X-interval(X).inf + 1; % input has to be positive
xs = [X.randPoint(100), X.randPoint(100, 'extreme')];
ys = nn.evaluate(xs);
[nn_red, Y_red] = nn.computeReducedNetwork(X, "InputCompression",true);

% check if points are contained
assert(all(zonotope(Y_red).contains(ys)));

% test completed
res = true;


end

% ------------------------------ END OF CODE ------------------------------
