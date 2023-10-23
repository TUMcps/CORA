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

resvec = [];

% test reading nnet network
nn = neuralNetwork.readNetwork('VertCAS_noResp_pra01_v9_20HU_200.nnet');
resvec(end+1) = true;

% test reading onnx network
nn = neuralNetwork.readNetwork('attitude_control_3_64_torch.onnx');
resvec(end+1) = true;

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
