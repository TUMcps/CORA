function P_out = convHull_(P,S,varargin)
% convHull_ - computes the convex hull of a polytope and another set or
%    point; to this end, we use the vertex representation of the polytope
%
% Syntax:
%    P_out = convHull_(P)
%    P_out = convHull_(P,S)
%
% Inputs:
%    P - polytope object
%    S - contSet object or point
%
% Outputs:
%    P_out - polytope enclosing the convex hull
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [4;3];
%    P = convHull(P1,P2);
%
%    figure; hold on
%    plot(P1,[1,2]);
%    plot(P2,[1,2]);
%    plot(P,[1,2],'k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull, conZonotope/convHull_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       26-November-2019
% Last update:   05-May-2020 (MW, standardized error message)
%                04-April-2022 (VK, adapted to polytope class)
%                27-July-2023 (MW, specialized method for 1D)
%                14-July-2024 (MW, support degenerate sets)
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% polytope is already convex (copy constructor)
if nargin == 1
    P_out = polytope(P);
    return; 
end

% ensure that numeric is second input argument
[P,S] = reorderNumeric(P,S);

% check dimensions
equalDimCheck(P,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < P.precedence
    P_out = convHull(S,P,varargin{:});
    return
end

% read out dimension
n = dim(P);

% convex hull with empty set is other set
if representsa_(P,'emptySet',eps) || representsa_(S,'emptySet',eps)
    P_out = polytope.empty(n);
    return
end

% convex hull with fullspace is fullspace
if representsa_(P,'fullspace',0)
    P_out = polytope.Inf(n);
    return
elseif representsa_(S,'fullspace',0)
    P_out = polytope.Inf(n);
    return
end


% special method for 1D case
if n == 1
    P_out = aux_convHull_1D(P,S);
    return
end

% numeric case: convex hull with a single vertex
if isnumeric(S)
    P_out = aux_convHull_vertices(P,S);
    return
end

% all other cases: try computing the vertex representation
try
    V2 = vertices(S);
catch ME
    throw(CORAerror('CORA:noops',P,S));
end
P_out = aux_convHull_vertices(P,V2);

end


% Auxiliary functions -----------------------------------------------------

function P_out = aux_convHull_vertices(P,V2)

% compute vertices of first polytope
V1 = vertices_(P,'lcon2vert');

% instantiate a polytope with the concatenated vertices
P_out = polytope([V1, V2]);

% reduce representation size
P_out = compact_(P_out,'V',1e-10);

% cannot be empty, otherwise premature exit
P_out.emptySet.val = false;

end

function P_out = aux_convHull_1D(P,S)

% compute vertices of both sets
V1 = vertices(P);
if isa(S,'contSet')
    V2 = vertices(S);
elseif isnumeric(S)
    V2 = S;
end

% check emptiness
if isempty(V1) || isempty(V2)
    P_out = polytope.empty(dim(P)); return
end

% compute convex hull of vertices
V = [min([V1,V2]),max([V1,V2])];
% check if there are actually two vertices
if withinTol(V(1),V(2))
    V = V(1);
end

% instantiate polytope
P_out = polytope(V);
P_out.emptySet.val = false;
P_out.fullDim.val = length(V) == 2;
P_out.bounded.val = ~any(isinf(V));

end

% ------------------------------ END OF CODE ------------------------------
