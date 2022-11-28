function res = testnn_requiredToolboxes
% testnn_requiredToolboxes - checks if required toolboxes are installed
% including the toolboxes for neural network verification
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

% Author:       Tobias Ladner
% Written:      02-November-2022
% Last update:  --
% Last revision:---

%------------- BEGIN CODE --------------

% test general toolboxes
res_partial(1) = test_requiredToolboxes;

% test nn verification specific toolboxes
res_partial(2) = logical(license('test', 'Neural_Network_Toolbox'));
if ~res_partial(2)
    disp('"Deep Learning Toolbox" missing!');
end

addons = matlab.addons.installedAddons;
res_partial(3) = any(strcmp("Deep Learning Toolbox Converter for ONNX Model Format", addons.Name));
if ~res_partial(3)
    disp('"Deep Learning Toolbox Converter for ONNX Model Format" missing!');
end

res = all(res_partial);

%------------- END OF CODE --------------
