function res = testnn_neuralNetwork_readONNXNetwork()
% testnn_neuralNetwork_readONNXNetwork - tests the readONNXNetwork function
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

% Authors:       Tobias Ladner
% Written:       28-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test reading basic network
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');

% test verbose output + input/output formats
nn = neuralNetwork.readONNXNetwork('controller_airplane.onnx', true, 'BC', 'BC');

% Reading network with custom layer
nn = neuralNetwork.readONNXNetwork('vnn_verivital_avgpool.onnx', false, 'BCSS');

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
