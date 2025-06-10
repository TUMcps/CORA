function res = test_nn_neuralNetwork_reduceGNNForNode()
% test_nn_neuralNetwork_reduceGNNForNode - tests the graph convolutional network 
%    functionality
%
% Syntax:
%    res = test_nn_neuralNetwork_reduceGNNForNode()
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
% Written:       26-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create graph (0-1, 0-3, 1-2, 1-4, 3-4)
adj_list = [; ...
    1, 1, 2, 2, 2, 3, 4, 4, 5, 5; ...
    2, 4, 1, 3, 5, 2, 1, 5, 2, 4; ...
    ];
G = graph(adj_list(1, :), adj_list(2, :), 'omitselfloops');
% add self loops
for i = 1:height(G.Nodes)
    G.addedge(i, i);
end

% create network
nn = neuralNetwork({ ...
    nnGNNLinearLayer([2, 3; 3, -4; 1, 2], [1; 2; 3]); ...
    nnGCNLayer(); ...
    nnReLULayer(); ...
    nnGNNLinearLayer([1, 2, 3; -4, 5, 2], [2; -1]); ...
    nnGCNLayer(); ...
    nnReLULayer(); ...
    nnGNNLinearLayer([1, 2; -4, -5; 2, 1], [-1; 4; -2]); ...
    });

% reduce for node
nn_red = nn.reduceGNNForNode(1,G);

% evaluate zonotope
c = [; ... % nodes are columns
    -2, 4, 1, -1, 3; ...
    5, 2, 3, 0, -1; ...
    ];
c = reshape(c, [], 1); % column vector
X = zonotope(c, eye(2*5)*0.01);

options = struct;
options.nn.graph = G;
nn.resetGNN();
Y = nn.evaluate(X, options);
nn.resetGNN();
Y_red = nn_red.evaluate(X, options);
assert(dim(Y_red) == 3)
assert(isequal(Y_red, project(Y,[1,6,11])))

% evaluate numeric
N = 500;
xs = randPoint(X, N);
nn.resetGNN();
ys = nn.evaluate(xs, options);
nn.resetGNN();
ys_red = nn_red.evaluate(xs, options);
assert(size(ys_red,1) == 3)
assert(compareMatrices(ys_red, ys([1,6,11],:)))

% check containment
res = all(Y.contains(ys, 'exact', 1e-12));

end

% ------------------------------ END OF CODE ------------------------------
