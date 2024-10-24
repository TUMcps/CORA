function [res] = testnn_nnVisualizationLayer()
% testnn_nnVisualizationLayer - tests constructor of nnVisualizationLayer
%
% Syntax:
%    res = testnn_nnVisualizationLayer()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       17-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create network
nn = neuralNetwork({ ...
    nnLinearLayer([2 -1; 4 2], [0.1; 0.2]); ...
    nnSigmoidLayer(); ...
    nnLinearLayer([5 -4; 1 -1; -1 2]); ...
    nnSigmoidLayer(); ...
});

% input set
X = zonotope([1;2], [1 0 1; 0 1 1]);
xs = X.randPoint(100);

% test propagation
ys = nn.evaluate(xs);
Y = nn.evaluate(X);

% add visualization
nnVis = nn.addVisualizationLayers();

% test propagation
ys = nnVis.evaluate(xs);
Y = nnVis.evaluate(X);

% close figures
close; close; close; close; close;

% test parameters
nnVis = nn.addVisualizationLayers([2 1], true);

% test propagation
ys = nnVis.evaluate(xs);
Y = nnVis.evaluate(X);

% close figures
close; close; close; close; close;

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
