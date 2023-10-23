function [res] = test_nn_neuralNetwork_readYMLNetwork()
% test_nn_neuralNetwork_readYMLNetwork - tests the readYMLNetwork function
%
% Syntax:
%    res = test_nn_neuralNetwork_readYMLNetwork()
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
% Written:       02-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Test 1
nn = neuralNetwork.readYMLNetwork('controllerACCtanh.yml');

res = true;

end

% ------------------------------ END OF CODE ------------------------------
