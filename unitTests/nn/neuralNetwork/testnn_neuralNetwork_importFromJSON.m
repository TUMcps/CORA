function res = testnn_neuralNetwork_importFromJSON()
% testnn_neuralNetwork_importFromJSON - unit test function for 
%     neuralNetwork/importFromJSON for ONNX network
%
% Syntax:
%    res = testnn_neuralNetwork_importFromJSON()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       10-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

nn = neuralNetwork.readONNXNetwork('attitude_control_3_64_torch.onnx');
nn = nn.getNormalForm();

nnJson = nn.exportAsJSON();
nn2 = neuralNetwork.importFromJSON(nnJson);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
