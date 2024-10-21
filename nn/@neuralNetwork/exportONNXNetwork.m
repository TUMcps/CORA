function res = exportONNXNetwork(nn, file_path, varargin)
% exportONNXNetwork - exports the given neural network to ONNX format
%
% Syntax:
%    res = exportONNXNetwork(nn,file_path)
%
% Inputs:
%    nn - neuralNetwork
%    file_path - path to file
%    varargin - additional name-value pairs for DLT ONNX export
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: exportONNXNetwork

% Authors:       Tobias Ladner
% Written:       30-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    file_path = './network.onnx';
end

% convert to DLT layers
nn_dlt = convertToDLToolboxNetwork(nn);

% export to ONNX
exportONNXNetwork(nn_dlt, file_path, varargin{:});

% result
res = true;

% ------------------------------ END OF CODE ------------------------------
