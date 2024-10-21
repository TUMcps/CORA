function res = test_nn_nnConv2DLayer()
% test_nn_nnConv2DLayer - tests constructor of nnConv2DLayer
%
% Syntax:
%    res = test_nn_nnConv2DLayer()
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
% See also: none

% Authors:       Tobias Ladner, Lukas Koller
% Written:       02-October-2023
% Last update:   22-February-2024 (LK, corrected false y_true)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple example
layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9]);

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],1);

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1]);

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1],[2 2],[2 2]);

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],1,[1 1 1 1],[2 2],[2 2],'testLayer');

% check evaluate
W = zeros(2,2,1,2);
W(:,:,:,1) = [1 -1; -1 2]; % filter 1
W(:,:,:,2) = [2 3; -1 -2]; % filter 2
b = [1; -2];
layer = nnConv2DLayer(W,b);
nn = neuralNetwork({layer});
n = 4;
nn.setInputSize([n,n,1]);

% check point
x = reshape(eye(n),[],1);
y = nn.evaluate(x);
y_true = [[...
    4 0 1
    0 4 0
    1 0 4],...
[
   -2 -3 -2
    1 -2 -3
   -2  1 -2]
];

assert(all(y == y_true(:)));

% check zonotope
X = zonotope(x,0.01 * eye(n*n));
Y = nn.evaluate(X);

assert(contains(Y,y));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
