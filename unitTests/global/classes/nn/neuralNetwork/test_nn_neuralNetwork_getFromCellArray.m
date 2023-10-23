function [res] = test_nn_neuralNetwork_getFromCellArray()
% test_nn_neuralNetwork_getFromCellArray - tests the getFromCellArray 
%    function supporting the legacy neuralNetworkOld constructor
%
% Syntax:
%    res = test_nn_neuralNetwork_getFromCellArray()
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
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

W = {rand(3, 2), rand(3, 3), rand(3, 2)};
b = {rand(3, 1), rand(3, 1), rand(3, 1)};
actFun = {'tanh', 'sigmoid', 'identity'};
nn = neuralNetwork.getFromCellArray(W, b, actFun);

W = {rand(3, 2), rand(3, 3), rand(3, 2)};
b = {rand(3, 1), rand(3, 1), rand(3, 1)};
actFun = 'relu';
nn = neuralNetwork.getFromCellArray(W, b, actFun);

try
    nn = neuralNetwork.getFromCellArray(rand(2, 2), rand(2, 1), actFun);
    % should have thrown an error
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
