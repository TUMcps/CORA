function [res] = test_nn_nnConv2DLayer()
% test_nn_nnConv2DLayer - tests constructor of 
%    nnConv2DLayer
%
% Syntax:
%    res = test_nn_nnConv2DLayer()
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
% Written:       02-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

try

% simple example
layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9]);
resvec(end+1) = true;

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],[1]);
resvec(end+1) = true;

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],[1],[1 1 1 1]);
resvec(end+1) = true;

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],[1],[1 1 1 1],[2 2],[2 2]);
resvec(end+1) = true;

layer = nnConv2DLayer([1 2 3; 4 5 6; 7 8 9],[1],[1 1 1 1],[2 2],[2 2],'testLayer');
resvec(end+1) = true;

% check evaluate
W = zeros(2,2,1,2);
W(:,:,:,1) = [1 -1; -1 2]; % filter 1
W(:,:,:,2) = [2 3; -1 -2]; % filter 2
b = [1 -2];
layer = nnConv2DLayer(W,b);
nn = neuralNetwork({layer});
n = 4;
nn.setInputSize([n,n,1]);

% check point
x = reshape(eye(n),[],1);
y = nn.evaluate(x);
y_true = [4 -3 1 -3 4 -3 1 -3 4 -2 4 -2 0 -2 4 -2 0 -2]';

resvec(end+1) = all(y == y_true);

% check zonotope
X = zonotope(x,0.01 * eye(n*n));
Y = nn.evaluate(X);

resvec(end+1) = contains(Y,y);

catch ME
    % failed test
    resvec(end+1) = false;
end


% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
