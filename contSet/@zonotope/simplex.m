function P = simplex(Z)
% simplex - enclose a zonotope by a simplex
%
% Syntax:
%    P = simplex(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    P - polytope object representing the simplex
%
% Example: 
%    Z = zonotope([1;0],[1 -1 0.5; 0 1 1]);
%    P = simplex(Z)
%    
%    figure; hold on; box on;
%    plot(Z);
%    plot(P, [1,2], 'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/interval, zonotope/polytope

% Authors:       Niklas Kochdumper
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct an n-dimensional standard simplex with origin 0
n = dim(Z);
V = eye(n+1);
B = gramSchmidt(ones(n+1,1));

P = polytope((B(:,2:end)'*V));

% scale the simplex so that it tightly encloses the zonotope    
A = P.A;
b = supremum(interval(A*Z));

P = polytope(A,b);

% ------------------------------ END OF CODE ------------------------------
