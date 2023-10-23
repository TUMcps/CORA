function V = vertices_(zB,varargin)
% vertices_ - Returns potential vertices of a zonotope bundle
%    WARNING: Do not use this function for high order zonotope bundles as
%    the computational complexity grows exponentially!
%
% Syntax:
%    V = vertices_(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    V - matrix storing the vertices 
%
% Example: 
%    Z1 = zonotope([1;1], [1 1; -1 1]);
%    Z2 = zonotope([-1;1], [1 0; 0 1]);
%    zB = zonoBundle({Z1,Z2});
%    V = vertices(zB);
% 
%    figure; hold on;
%    plot(zB);
%    plot(V(1,:),V(2,:),'k.','MarkerSize',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices, polytope

% Authors:       Matthias Althoff
% Written:       18-August-2016 
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename vertices_)

% ------------------------------ BEGIN CODE -------------------------------

%obtain polytope
P = polytope(zB);

%obtain vertices (input check in polytope-function)
V = vertices(P);

% return correct dimension
if isempty(V)
    V = double.empty(dim(zB),0);
end

% ------------------------------ END OF CODE ------------------------------
