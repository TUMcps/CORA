function [res] = testnn_neuralNetwork_readONNXNetwork()
% testnn_neuralNetwork_readONNXNetwork - tests the readONNXNetwork funciton
%
% Syntax:  
%    res = testnn_neuralNetwork_readONNXNetwork()
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

% Author:       Tobias Ladner
% Written:      28-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% Test 1
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');

% Test 2
nn = neuralNetwork.readONNXNetwork('controller_airplane.onnx', true, 'BC', 'BC');

end

%------------- END OF CODE --------------

