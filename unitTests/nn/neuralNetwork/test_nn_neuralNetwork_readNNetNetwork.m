function [res] = test_nn_neuralNetwork_readNNetNetwork()
% test_nn_neuralNetwork_readNNetNetwork - tests the readNNetNetwork 
%    function
%
% Syntax:
%    res = test_nn_neuralNetwork_readNNetNetwork()
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
nn = neuralNetwork.readNNetNetwork('VertCAS_noResp_pra01_v9_20HU_200.nnet');

res = true;

end

% ------------------------------ END OF CODE ------------------------------
