function res = isequal(pHA1,pHA2,varargin)
% isequal - checks if two parallel hybrid automata are equal
%
% Syntax:
%    res = isequal(pHA1,pHA2)
%    res = isequal(pHA1,pHA2,tol)
%
% Inputs:
%    pHA1 - parallelHybridAutomaton object
%    pHA2 - parallelHybridAutomaton object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      10-January-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{pHA1,'att','parallelHybridAutomaton'};
                {pHA2,'att','parallelHybridAutomaton'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check number of components
if length(pHA1.components) ~= length(pHA2.components)
    res = false; return
end

% loop over components
for i=1:length(pHA1.components)
    
    % compare hybrid automata
    if ~isequal(pHA1.components{i},pHA2.components{i},tol)
        res = false; return
    end

    % compare number input binds
    if any(size(pHA1.bindsInputs{i}) ~= size(pHA2.bindsInputs{i}))
        res = false; return
    end

    % compare input binds
    if any(any(pHA1.bindsInputs{i} ~= pHA2.bindsInputs{i}))
        res = false; return
    end

end

% all checks ok
res = true;

%------------- END OF CODE --------------