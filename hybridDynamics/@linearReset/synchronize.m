function linReset_sync = synchronize(linResets,varargin)
% synchronize - synchronize linear reset functions of equal pre-state and
%    post-state dimensions by adding and concatenating their respective
%    matrices, which corresponds to the reset functions being evaluated
%    synchronously
%
% Syntax:
%    linReset_sync = synchronize(linResets)
%    linReset_sync = synchronize(linResets,idStates)
%
% Inputs:
%    linResets - class array of linearReset objects
%    idStates - indices of states that are mapped by identity
%
% Outputs:
%    linReset_sync - linearReset object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/synchronize

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       04-April-2022
% Last update:   01-July-2022
%                14-January-2023 (MW, handle states unaffected by sync)
% Last revision: 10-October-2024 (MW, moved from transition class)

% ------------------------------ BEGIN CODE -------------------------------

% set default value for identity mapping
preStateDim = linResets(1).preStateDim;
postStateDim = linResets(1).postStateDim;
idStates = setDefaultValues({[]},varargin);

% all reset functions need to have the same mapping
preStateDims = arrayfun(@(r) r.preStateDim,linResets,'UniformOutput',true);
postStateDims = arrayfun(@(r) r.postStateDim,linResets,'UniformOutput',true);
if any(postStateDims ~= postStateDims(1)) || any(preStateDims ~= preStateDims(1))
    throw(CORAerror('CORA:wrongValue','first',...
        'All reset functions must map the same spaces: R^n -> R^m.'));
end
inputDim = linResets(1).inputDim;

% initialize state/input matrix and constant offset
A = zeros(postStateDim,preStateDim);
B = zeros(postStateDim,inputDim);
c = zeros(postStateDim,1);

% since all resets have been projected to higher dimensions before,
% we can just add them here
for i=1:length(linResets)
    A = A + linResets(i).A;
    B = B + linResets(i).B;
    c = c + linResets(i).c;
end

% insert identity mapping
id = zeros(postStateDim,1);
id(idStates) = 1;
A = A + diag(id);

% init synchronized reset
linReset_sync = linearReset(A,B,c);

% ------------------------------ END OF CODE ------------------------------
