function res = testnn_neuralNetwork_addVisualizationLayers()
% testnn_neuralNetwork_addVisualizationLayers - unit test function for 
%     neuralNetwork/addVisualizationLayers
%
% Syntax:
%    res = testnn_neuralNetwork_addVisualizationLayers()
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
% See also: neuralNetwork/addVisualizationLayers

% Authors:       Tobias Ladner
% Written:       02-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

nn = neuralNetwork({ ...
    nnLinearLayer([-1 2; 1 3],[1;2]); ...
    nnSigmoidLayer();
    nnLinearLayer([-1 2; -1 3],[1;-2]); ...
    nnLeakyReLULayer();
});

nn_vis = nn.addVisualizationLayers();

% test evaluate

% test point
X = zonotope([1;2],[0.1 -0.02; 0.03 0.2]);
xs = X.randPoint(20);
ys = nn.evaluate(xs);
ys_vis = nn_vis.evaluate(xs);
assert(isequal(ys,ys_vis));

% test if figures were created
fids = 10001:10005;
assert(all(ishandle(fids)));

% test zonotope
Y = nn.evaluate(X);
Y_vis = nn_vis.evaluate(X);
assert(isequal(Y,Y_vis));

% test if figures were created
assert(all(ishandle(fids)));

% test interval
X = interval(X);
Y = nn.evaluate(X);
Y_vis = nn_vis.evaluate(X);
assert(isequal(Y,Y_vis));

% test if figures were created
assert(all(ishandle(fids)));

% test if all figures have three plots (numeric, zonotope, interval)
for fid=fids
    figure(fid)
    assertLoop(length(gca().Children) == 3,fid);
end

% close figures
close(fids)

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
