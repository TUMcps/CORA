function [V,empty] = priv_vertices_1D(A,b,Ae,be)
% priv_vertices_1D - computes the vertices of a 1D polytope
%    note: unbounded 1D polytopes are supported
%
% Syntax:
%    [V,empty] = priv_vertices_1D(A,b,Ae,be,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%
% Outputs:
%    V - vertex representation
%    empty - true/false whether polytope is empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set tolerance
tol = 1e-12;

% compute minimal representation
[A,b,Ae,be,empty] = priv_compact_1D(A,b,Ae,be,tol);

% if there is no point, P is already empty
if empty
    V = [];
    return
end

% normalize constraints
[A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'A');

if ~isempty(A)

    % check boundedness from below
    Aisminus1 = withinTol(A,-1,tol);
    if any(Aisminus1)
        % bounded from below
        V = -b(Aisminus1);
    else
        % unbounded toward -Inf
        V = -Inf;
    end

    % check boundedness from above
    Ais1 = withinTol(A,1,tol);
    if any(Ais1)
        % bounded from above (add only if not a duplicate)
        if ~withinTol(V,b(Ais1))
            V = [V, b(Ais1)];
        end
    else
        % unbounded toward +Inf
        V = [V, Inf];
    end

elseif ~isempty(Ae)
    % due to minHRep call above, we should only have one equality here
    V = be / Ae;

else
    throw(CORAerror('CORA:specialError',...
        'Error in vertex computation of 1D polytope.'));

end

% ------------------------------ END OF CODE ------------------------------
