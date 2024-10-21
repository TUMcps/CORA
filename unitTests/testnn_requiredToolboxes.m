function res = testnn_requiredToolboxes
% testnn_requiredToolboxes - checks if required toolboxes are installed
%    including the toolboxes for neural network verification.
%    Unit tests requiring these additional toolboxes are indicated by
%    'testnn' instead of 'test_nn.
%
% Syntax:
%    res = testnn_requiredToolboxes
%
% Inputs:
%    nolice
%
% Outputs:
%    res - boolean 
%
% Example: 
%

% Authors:       Tobias Ladner
% Written:       02-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test general toolboxes
test_requiredToolboxes;

% test nn verification specific toolboxes

% Deep Learning Toolbox
assert(logical(license('test', 'Neural_Network_Toolbox')),'Deep Learning Toolbox missing!');

% ONNX support package
text = "Deep Learning Toolbox Converter for ONNX Model Format";
addons = matlabshared.supportpkg.getInstalled;
assert(~isempty(addons) && ismember(text, {addons.Name}),'Deep Learning Toolbox Converter for ONNX Model Format is missing!');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
