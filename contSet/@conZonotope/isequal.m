function res = isequal(cZ,S,varargin)
% isequal - checks if a constrained zonotope represents the same set as
%    another set
%
% Syntax:
%    res = isequal(cZ,S)
%    res = isequal(cZ,S,tol)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    c = [0;0]; G = [1 0 1;0 1 1];
%    A = [1 1 1]; b = 4;
%    cZ = conZonotope([c,G],A,b);
%
%    isequal(cZ,cZ + [0;eps],0)
%    isequal(cZ,cZ + [0;eps],1e-14)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-December-2022
% Last update:   16-December-2023 (MW, support comparison to polytope)
%                18-July-2024 (MW, centers may be different)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% set default tolerance
tol = setDefaultValues({1e-12},varargin);

% check input arguments
inputArgsCheck({ {cZ,'att',{'conZonotope','numeric'}}, ...
                 {S,'att',{'contSet','numeric'}}, ...
                 {tol,'att','numeric',{'nonempty','scalar','finite','nonnegative'}}});

% ensure that numeric is second input argument
[cZ,S] = reorderNumeric(cZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cZ.precedence
    res = isequal(S,cZ,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(cZ,S,true)
    res = false;
    return
end

% check empty cases
cZ_empty = representsa_(cZ,'emptySet',tol);
S_empty = representsa_(S,'emptySet',tol);
if cZ_empty && S_empty
    res = true;
    return
elseif xor(cZ_empty,S_empty)
    res = false;
    return
end

% two constrained zonotopes
if isa(S,'conZonotope')
    res = aux_isequal_conZonotope(cZ,S,tol);
    return
end

% second set: zonotope or interval
if isa(S,'zonotope') || isa(S,'interval')
    % first set: conZonotope without constraints
    if representsa_(cZ,'zonotope',tol)
        res = isequal(zonotope(cZ),zonotope(S),tol);
    else
        res = false;
    end
    return
end

throw(CORAerror('CORA:noops',cZ,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_conZonotope(cZ,S,tol)

% check point cases
[cZ_isPoint, cZ_p] = representsa_(cZ,'point',tol);
[S_isPoint, S_p] = representsa_(S,'point',tol);
if cZ_isPoint && S_isPoint
    res = all(withinTol(cZ_p, S_p, tol));
    return
elseif (~cZ_isPoint && S_isPoint) || ...
        (cZ_isPoint && ~S_isPoint)
    res = false;
    return
end

% check if both are represented completely equally
if all(withinTol(cZ.c,S.c,tol)) ...
        && compareMatrices([cZ.G; cZ.A],[S.G; S.A],tol) ...
        && compareMatrices([cZ.A cZ.b],[S.A S.b],tol)
    res = true;
    return
end

% general method: check vertices (computationally expensive!)
res = compareMatrices(vertices(cZ),vertices(S),tol);

end

% ------------------------------ END OF CODE ------------------------------
