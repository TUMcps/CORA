function [res] = test_nn_neuralNetwork_readSherlockNetwork()
% test_nn_neuralNetwork_readSherlockNetwork - tests the readSherlockNetwork
%    function
%
% Syntax:
%    res = test_nn_neuralNetwork_readSherlockNetwork()
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
nn = neuralNetwork.readSherlockNetwork('controllerAirplane.sherlock');

res = true;

end

% ------------------------------ END OF CODE ------------------------------
