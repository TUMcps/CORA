function res = isequal(P,S,varargin)
% isequal - check ifs two polytopes are equal
%
% Syntax:
%    res = isequal(P,S)
%    res = isequal(P,S,tol)
%
% Inputs:
%    P - polytope object 
%    S - contSet object 
%    tol - (optional) tolerance
%
% Outputs:
%    res - result of comparison
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([0 -1; 0 1;-1 0; 1 0; -1 -1],[2;3;2;3;2]);
%    res = isequal(P1,P2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-April-2023 (MW, migrated code from eq)
% Last update:   27-July-2023 (MW, handle 1D case)
%                16-December-2023 (MW, comparison to numeric)
%                02-December-2023 (MW, fix fully empty polytopes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
tol = setDefaultValues(1e-9,varargin);

% check dimensions
try
    equalDimCheck(P,S);
catch ME
    % P and S have a different ambient dimension
    res = false; return
end

% handle fully empty polytope objects
if representsa_(P,'fullspace',0)
    res = representsa_(S,'fullspace',0); return;
elseif representsa_(S,'fullspace',0)
    % if P were fullspace, we would have entered to if-branch above
    res = false; return
end

% one-dimensional case
if dim(P) == 1
    res = aux_1D(P,S,tol); return
end

% quick check for polytopes: check emptiness, boundedness, and degeneracy

% different set representations
if isa(S,'polytope') && aux_diffProperties(P,S)
    res = false; return
end

if isnumeric(S)
    % avoid comparison with matrices
    if size(S,2) > 1
        throw(CORAerror('CORA:noops',P,S));
    end

    % S is a single point
    if ~representsa_(P,'point',tol)
        res = false;
    else
        % compute vertex
        V = vertices(P);
        res = all(withinTol(V,S,tol));
    end

else 

    % check if P and S contain each other
    res = contains_(P,S,'exact',tol) && contains_(S,P,'exact',tol);
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_1D(P,S,tol)

    % compute vertices (fast)
    V1 = vertices_(P,'lcon2vert');
    if isnumeric(S)
        V2 = S;
    else
        V2 = vertices(S);
    end

    % compare
    if xor(any(isinf(V1)),any(isinf(V2)))
        res = false;
    elseif size(V1,2) ~= size(V2,2) 
        res = false;
    elseif any(isinf(V1))
        res = all(V1 == V2);
    else
        res = compareMatrices(V1,V2,tol);
    end

end

function res = aux_diffProperties(P1,P2)
% returns true if two polytopes P1 and P2 have the different properties (if
% given), e.g., P1 is bounded but P2 is unbounded
% if properties unknown, we return false

% assume same properties (or unknown)
res = false;

% boundedness
if ~isempty(P1.bounded.val) && ~isempty(P2.bounded.val) ...
        && xor(P1.bounded.val,P2.bounded.val)
    res = true; return
end

% emptiness
if ~isempty(P1.emptySet.val) && ~isempty(P2.emptySet.val) ...
        && xor(P1.emptySet.val,P2.emptySet.val)
    res = true; return
end

% degeneracy
if ~isempty(P1.fullDim.val) && ~isempty(P2.fullDim.val) ...
        && xor(P1.fullDim.val,P2.fullDim.val)
    res = true; return
end

end

% ------------------------------ END OF CODE ------------------------------
