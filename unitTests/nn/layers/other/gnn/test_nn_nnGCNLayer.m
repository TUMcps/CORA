function res = test_nn_nnGCNLayer()
% test_nn_nnGCNLayer - tests the graph convolutional network functionality
%
% Syntax:
%    res = test_nn_nnGCNLayer()
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
    nnGCNLayer(); ...
});

% evaluate zonotope
c = [; ... % nodes are columns
    -2, 4, 1, -1, 3; ...
    5, 2, 3, 0, -1; ...
    ]';
c = reshape(c, [], 1); % column vector
X = zonotope(c, eye(2*5)*0.01);

options = struct;
options.nn.graph = G;
nn.resetGNN();
Y = nn.evaluate(X, options);

% evaluate numeric
N = 500;
xs = randPoint(X, N);
ys = nn.evaluate(xs, options);

% check containment
assert(all(Y.contains(ys, 'exact', 1e-12)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
