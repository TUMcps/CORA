function res = isequal(Z,S,varargin)
% isequal - checks if a zonotope is equal to another set or point
%    note: no deletion of aligned generators since this is quite costly
%
% Syntax:
%    res = isequal(Z,S)
%    res = isequal(Z,S,tol)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object, numeric
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0 2; 2 0 1]);
%    Z2 = zonotope(zeros(2,1),[1 2 0; 2 1 0]);
%    isequal(Z1,Z2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   09-June-2020 (include re-ordering of generators)
%                13-November-2022 (MW, integrate modular check)
%                05-October-2024 (MW, fix 1D case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{Z,'att',{'zonotope','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% ensure that numeric is second input argument
[Z,S] = reorderNumeric(Z,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < Z.precedence
    res = isequal(S,Z,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(Z,S,true)
    res = false;
    return
end

if isa(S,'zonotope')
    res = aux_isequal_zonotope(Z,S,tol);
    return
end

% other sets: convert to zonotope (exact)
if isa(S,'interval') || isnumeric(S)
    S = zonotope(S);
    res = aux_isequal_zonotope(Z,S,tol);
    return
end

throw(CORAerror('CORA:noops',Z,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_zonotope(Z,S,tol)

% init result
res = false;

% compare centers (quick check)
if ~all(withinTol(Z.c,S.c,tol))
    return
end

% obtain minimal representation
G1 = compact_(Z,'all',tol).G;
G2 = compact_(S,'all',tol).G;

% compare number of generators
if size(G1,2) ~= size(G2,2)
    return
end

% compare generator matrices: must match with no remainder, order is
% irrelevant, sign may be inverted
res = compareMatrices(G1,G2,tol,'equal',false,false);

end

% ------------------------------ END OF CODE ------------------------------
