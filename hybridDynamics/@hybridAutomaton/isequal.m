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

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   21-May-2023 (MW, extend to arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% check array length
if any(size(HA1) ~= size(HA2))
    res = false; return
end

% loop over all individual transition objects
[r,c] = size(HA1);
res = false(r,c);
for i=1:r
    for j=1:c
        res(i,j) = aux_isequal(HA1(i,j),HA2(i,j),tol);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal(HA1,HA2,tol)

% check number of locations
if length(HA1.location) ~= length(HA2.location)
    res = false; return
end

% loop over locations
for i=1:length(HA1.location)
    
    % compare locations
    if ~isequal(HA1.location(i),HA2.location(i),tol)
        res = false; return
    end

end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
