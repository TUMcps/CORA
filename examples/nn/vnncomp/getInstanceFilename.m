function [instanceFilename,modelName,vnnlibName] = ...
    getInstanceFilename(benchName,modelPath,vnnlibPath)
% getInstanceFilename - create a filename for storing the parsed neural 
%   network and vnnlib specification as a .mat file.
%
% Syntax:
%    [instanceFilename,modelName,vnnlibName] = ...
%       getInstanceFilename(benchName,modelPath,vnnlibPath)
%
% Inputs:
%    benchName - name of the benchmark
%    modelPath - path to the .onnx-file
%    vnnlibPath - path to the .vnnlib-file
%
% Outputs:
%    instanceFilename - filename (unique for this instance)
%    modelName - name of the .onnx-file
%    vnnlibName - name of the .vnnlib-file
%
% References:
%    [1] VNN-COMP'24
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

modelName = regexp(modelPath,'([^/]+)(?=\.onnx$)','match');
vnnlibName = regexp(vnnlibPath,'([^/]+)(?=\.vnnlib$)','match');
instanceFilename = [benchName '_' modelName{1} '_' vnnlibName{1} '.mat'];

end

% ------------------------------ END OF CODE ------------------------------
