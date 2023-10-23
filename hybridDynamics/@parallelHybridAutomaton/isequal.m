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
inputArgsCheck({{pHA1,'att','parallelHybridAutomaton'};
                {pHA2,'att','parallelHybridAutomaton'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check array length
if any(size(pHA1) ~= size(pHA2))
    res = false; return
end

% loop over all individual transition objects
[r,c] = size(pHA1);
res = false(r,c);
for i=1:r
    for j=1:c
        res(i,j) = aux_isequal(pHA1(i,j),pHA2(i,j),tol);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal(pHA1,pHA2,tol)

% check number of components
if length(pHA1.components) ~= length(pHA2.components)
    res = false; return
end

% loop over components
for i=1:length(pHA1.components)

    % check for emptiness
    if xor(isemptyobject(pHA1),isemptyobject(pHA2))
        res = false; return
    end
    
    if ~isemptyobject(pHA1)
        % compare hybrid automata
        if ~isequal(pHA1.components(i),pHA2.components(i),tol)
            res = false; return
        end
    
        % compare number of input binds
        if any(size(pHA1.bindsInputs{i}) ~= size(pHA2.bindsInputs{i}))
            res = false; return
        end
    
        % compare input binds
        if any(any(pHA1.bindsInputs{i} ~= pHA2.bindsInputs{i}))
            res = false; return
        end
    end

end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
