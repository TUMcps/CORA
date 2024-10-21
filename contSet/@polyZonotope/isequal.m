function res = isequal(pZ,S,varargin)
% isequal - checks if a polynomial zonotope is equal to another set or
%    point
%
% Syntax:
%    res = isequal(pZ,S)
%    res = isequal(pZ,S,tol)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object, numeric
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([0;0],[1 0 1;0 -1 1],[0.4 0;0.1 1],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 1 0;1 0 -1],[0 0.4;1 0.1],[2 1 0;1 0 1]);
%    isequal(pZ1,pZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isequal

% Authors:       Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{pZ,'att',{'polyZonotope','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% ensure that numeric is second input argument
[pZ,S] = reorderNumeric(pZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pZ.precedence
    res = isequal(S,pZ,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(pZ,S,true)
    res = false;
    return
end

if isa(S,'polyZonotope')
    res = aux_isequal_polyZonotope(pZ,S,tol);
    return
end

% for zonotopes and interval, use exact conversion (fast)
if isa(S,'zonotope') || isa(S,'interval')
    S = polyZonotope(S);
    res = aux_isequal_polyZonotope(pZ,S,tol);
    return
end

throw(CORAerror('CORA:noops',pZ,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_polyZonotope(pZ,S,tol)

% assume false
res = false;

% remove redundancies in representation
pZ = compact_(pZ,'all',eps);
S = compact_(S,'all',eps);

% compare number of generators (quick check)
if size(pZ.G,2) ~= size(S.G,2) || size(pZ.GI,2) ~= size(S.GI,2)
    return 
end

% compare identifier vectors
temp1 = sort(pZ.id); temp2 = sort(unique([pZ.id;S.id]));
if length(temp1) ~= length(temp2) || ~all(temp1 == temp2)
    return;
elseif ~all(pZ.id == S.id)
    [~,E1,E2] = mergeExpMatrix(pZ.id,S.id,pZ.E,S.E);
else
    E1 = pZ.E; E2 = S.E;
end

% jointly compare dependent generators and exponent matrices
if ~compareMatrices([pZ.G;E1],[S.G;E2], tol)
    return
end
if ~compareMatrices(pZ.GI,S.GI, tol, 'equal', false, false)
    return
end
if ~all(withinTol(pZ.c, S.c, tol))
    return
end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
