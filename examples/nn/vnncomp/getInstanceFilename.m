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
% Last update:   02-June-2026 (BK, multi-network path support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% match the .onnx stem anywhere in the string, excluding Python-list chars so
% multi-network paths like "[('f','onnx/foo.onnx'),('g',...)]" also work
modelName = regexp(modelPath,'([^/''"() ]+)(?=\.onnx)','match');
vnnlibName = regexp(vnnlibPath,'([^/]+)(?=\.vnnlib)','match');
instanceFilename = [benchName '_' modelName{1} '_' vnnlibName{1} '.mat'];

end

% ------------------------------ END OF CODE ------------------------------
