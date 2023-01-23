function res = isequal(HA1,HA2,varargin)
% isequal - checks if two hybrid automata are equal
%
% Syntax:  
%    res = isequal(HA1,HA2)
%    res = isequal(HA1,HA2,tol)
%
% Inputs:
%    HA1 - hybridAutomaton object
%    HA2 - hybridAutomaton object
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
inputArgsCheck({{HA1,'att','hybridAutomaton'};
                {HA2,'att','hybridAutomaton'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check number of locations
if length(HA1.location) ~= length(HA2.location)
    res = false; return
end

% loop over locations
for i=1:length(HA1.location)
    
    % compare locations
    if ~isequal(HA1.location{i},HA2.location{i},tol)
        res = false; return
    end

end

% all checks ok
res = true;

%------------- END OF CODE --------------