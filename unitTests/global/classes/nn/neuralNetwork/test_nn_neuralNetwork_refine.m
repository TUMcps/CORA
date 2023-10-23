function res = test_nn_neuralNetwork_refine()
% test_nn_neuralNetwork_refine - unit test function for 
%     neuralNetwork/evaluate: sanity check of points
%
% Syntax:
%    res = test_nn_neuralNetwork_refine()
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

% Authors:       Tobias Ladner
% Written:       13-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create network
nn = neuralNetwork({ ...
    nnLinearLayer([2 3; 4 5]); ...
    nnSigmoidLayer(); ...
    nnLinearLayer([-1 5; 2 -3]); ...
    nnSigmoidLayer(); ...
});

% evaluate point
x = [5;2];
y = nn.evaluate(x);

% add uncertainty
X = polyZonotope(x, eye(2)*0.01);

% test refine functions ---

nn.reset();
nn.evaluate(X);
nn.refine()
resvec(end+1) = true;

nn.reset();
nn.evaluate(X);
nn.refine(2)
resvec(end+1) = true;

nn.reset();
nn.evaluate(X);
nn.refine(5, "layer")
nn.refine(5, "neuron")
nn.refine(5, "all")
resvec(end+1) = true;

nn.reset();
nn.evaluate(X);
nn.refine(5, "layer", "approx_error");
nn.refine(5, "layer", "sensitivity", x);
nn.refine(5, "layer", "both", x);
nn.refine(5, "layer", "random");
resvec(end+1) = true;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
