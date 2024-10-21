function res = isequal(C,S,varargin)
% isequal - checks if two capsules are equal
%
% Syntax:
%    res = isequal(C,S)
%    res = isequal(C,S,tol)
%
% Inputs:
%    C - capsule object
%    S - capsule object, numeric
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    C1 = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    C2 = capsule([1; 0; 0], [0.5; -1; 1], 0.5);
%    res = isequal(C1,C2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   12-March-2021 (MW, add dimension mismatch)
%                01-June-2022 (MW, correct dimension mismatch case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{C,'att',{'capsule','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% ensure that numeric is second input argument
[C,S] = reorderNumeric(C,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < C.precedence
    res = isequal(S,C,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(C,S,true)
    res = false;
    return
end

% capsule-capsule case
if isa(S,'capsule')
    res = aux_isequal_capsule(C,S,tol);
    return
end

throw(CORAerror('CORA:noops',C,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_capsule(C,S,tol)

% check emptiness
C1_empty = representsa_(C,'emptySet',tol);
C2_empty = representsa_(S,'emptySet',tol);
if xor(C1_empty,C2_empty)
    res = false; return
elseif C1_empty && C2_empty
    res = true; return
end

% check for equality
res = all(abs(center(C) - center(S)) < tol) && ... % center
    all(abs(C.g - S.g) < tol) && ... % generator
    abs(C.r - S.r) < tol; % radius

end

% ------------------------------ END OF CODE ------------------------------
