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

resvec = [];

% test reading basic network
nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');
resvec(end+1) = true;

% test verbose output + input/output formats
nn = neuralNetwork.readONNXNetwork('controller_airplane.onnx', true, 'BC', 'BC');
resvec(end+1) = true;

% Reading network with custom layer
nn = neuralNetwork.readONNXNetwork('vnn_verivital_avgpool.onnx', false, 'BCSS');
resvec(end+1) = true;

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
