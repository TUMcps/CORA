function P = convHull(P,S)
% convHull - computes the convex hull of two polytopes
%
% Syntax:  
%    P = convHull(P1,S)
%
% Inputs:
%    P - mptPolytope object
%    S - contSet object
%
% Outputs:
%    P - mptPolytope enclosing the convex hull of P and S
%
% Example: 
%    P1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [4;3];
%
%    P = convHull(P1,P2);
%
%    figure; hold on;
%    plot(P1,[1,2],'FaceColor'.'r');
%    plot(P2,[1,2],'FaceColor','b');
%    plot(P,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:        Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 1
    return;
end

% find a polytope object
[P,S] = findClassArg(P,S,'mptPolytope');

% different cases depending on the class of the summand
if isnumeric(S)
    
    V2 = S;

elseif isa(S,'mptPolytope') || isa(S,'interval') || ...
       isa(S,'conZonotope') || isa(S,'zonoBundle') || ...
       isa(S,'zonotope')

    % compute vertices
    V2 = vertices(S);

else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',P,S));
end

% comptue vertices of first polytope
V1 = vertices(P);

% compute convex hull
V = [V1,V2];
K = convhulln(V');
indices = unique(K);

% instantiate resulting polytope
P = mptPolytope(V(:,indices)');

%------------- END OF CODE --------------