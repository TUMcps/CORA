function res = isequal(zB,S,varargin)
% isequal - checks if a zonotope bundle is equal to another set or point
%
% Syntax:
%    res = isequal(zB,S)
%    res = isequal(zB,S,tol)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object, numeric
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    Z1 = zonotope([1;1], [1 1; -1 1]);
%    Z2 = zonotope([-1;1], [1 0; 0 1]);
%    zB = zonoBundle({Z1,Z2});
%    isequal(zB,zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       17-September-2019
% Last update:   08-July-2024 (AK, implemented exact checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{zB,'att',{'zonoBundle','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% ensure that numeric is second input argument
[zB,S] = reorderNumeric(zB,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < zB.precedence
    res = isequal(S,zB,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(zB,S,true)
    res = false;
    return
end

if isa(S,'zonoBundle')
    res = aux_isequal_zonoBundle(zB,S,tol);
    return
end

% general idea: check using mutual containment
if isa(S,'contSet')
    res = zB.contains_(S,'exact',tol,0,false,false) && S.contains_(zB,'exact',tol,0,false,false);
    return
end

throw(CORAerror('CORA:noops',zB,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_zonoBundle(zB,S,tol)

res = true;
% quick check: all individual zonotopes are equal
if zB.parallelSets == S.parallelSets
    for z=1:zB.parallelSets
        if ~isequal(zB.Z{z},S.Z{z})
            res = false;
            break
        end
    end
end
if res
    return
end

% general check using mutual containment
res = zB.contains(S, 'exact', tol) && S.contains(zB, 'exact', tol);

end

% ------------------------------ END OF CODE ------------------------------
