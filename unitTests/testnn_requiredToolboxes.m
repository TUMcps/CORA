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
res_partial(1) = test_requiredToolboxes;

% test nn verification specific toolboxes

% Deep Learning Toolbox
res_partial(2) = logical(license('test', 'Neural_Network_Toolbox'));
if ~res_partial(2)
    disp('"Deep Learning Toolbox" missing!');
end

% ONNX support package
text = "Deep Learning Toolbox Converter for ONNX Model Format";
addons = matlabshared.supportpkg.getInstalled;
res_partial(3) = ~isempty(addons) && ismember(text, {addons.Name});
if ~res_partial(3)
    disp('"Deep Learning Toolbox Converter for ONNX Model Format" missing!');
end

res = all(res_partial);

% ------------------------------ END OF CODE ------------------------------
