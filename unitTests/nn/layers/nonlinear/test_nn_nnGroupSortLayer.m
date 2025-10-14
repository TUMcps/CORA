function res = test_nn_nnGroupSortLayer()
% test_nn_nnGroupSortLayer - tests the ReLU layer
%
% Syntax:
%    res = test_nn_nnGroupSortLayer()
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

% Authors:       Lukas Koller
% Written:       18-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Specify a group size.
groupSize = 3;

% Initialize a GroupSort Layer.
layer = nnGroupSortLayer(groupSize);

% Test attribute.
assert(layer.groupSize == groupSize);

% Specify number of input dimensions.
n = 11;
% Specify a batch size.
bSz = 5;
% Generate a random input.
x = rand([n bSz]);

% Enable backpropagation, to check the stored permutation.
options.nn.train.backprop = true;

% Apply groupSort.
y = layer.evaluate(x,options);

% Check the dimensionality of the output.
assert(all(size(y) == size(x)));

% Check if all dimensions are sorted.
for i=1:bSz % Iterate over the batch elements.
    for j=1:groupSize:n % Iterate over the groups.
        % Obtain the current group.
        gji = y(j:min(j+groupSize-1,end),i);
        [~,idx] = sort(gji,1,'ascend');
        assert(all(idx' == 1:length(gji)));
    end
end

% Obtain the stored permutation.
permIdx = layer.backprop.store.permIdx;
% Check if permutation is correct.
assert(all(y == x(permIdx),'all'));

res = true;


% ------------------------------ END OF CODE ------------------------------
