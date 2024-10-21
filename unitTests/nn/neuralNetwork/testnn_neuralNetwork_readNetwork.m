function res = testnn_neuralNetwork_readNetwork()
% testnn_neuralNetwork_readNetwork - tests the readNetwork function
%
% Syntax:
%    res = testnn_neuralNetwork_readNetwork()
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

% test reading nnet network
nn = neuralNetwork.readNetwork('VertCAS_noResp_pra01_v9_20HU_200.nnet');

% test reading onnx network
nn = neuralNetwork.readNetwork('attitude_control_3_64_torch.onnx');

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
