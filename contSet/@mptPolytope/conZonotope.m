function res = conZonotope(obj)
% conZonotope - convert a mptPolytope object into a constrained zonotope
%               object
%
% Syntax:  
%    res = conZonotope(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - c-zonotope object
%
% Example: 
%    A = [-1 0;0 -1;1 1];
%    b = [1;1;1];
%    poly = mptPolytope(A,b);
%    cZono = conZonotope(poly);
%
%    hold on
%    plotFilled(cZono,[1,2],'r','EdgeColor','none')
%    plot(poly,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% Implementation according to Theorem 1 form [1]:

% extract object properties
A = obj.P.A;
b = obj.P.b;

% calculate the vertices of the polytope
v = vertices(obj);
V = get(v,'V');

% calculate a bounding box for the constrained zonotope and convert it to a
% zonotope
bb = interval(min(V,[],2),max(V,[],2));
zono = zonotope(bb);

Z = get(zono,'Z');
c = Z(:,1);
G = Z(:,2:end);

% Calculate the lower bound sigma for A*x \in [sigma,b] (Theorem 1 in [1])
sigma = min(A*V,[],2);

% Construct constrained zonotope object according to equation (21) in [1]
G_ = [G, zeros(size(G,1),size(A,1))];
A_ = [A*G, diag((sigma-b)./2)];
b_ = (b+sigma)./2 - A*c;

res = conZonotope([c,G_],A_,b_);


%------------- END OF CODE --------------