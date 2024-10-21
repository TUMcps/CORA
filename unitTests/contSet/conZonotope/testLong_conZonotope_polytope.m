function res = testLong_conZonotope_polytope
% testLong_conZonotope_polytope - unit test function for
%    conversion between constrained zonotopes and polytopes
%
% Syntax:
%    res = testLong_conZonotope_polytope
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 10;

for j=1:nrTests

    % random dimension
    n = randi([2,3]);

    % Generate random polytope vertices
    points = rand(n,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a polytope object from the vertices
    P = polytope(V);

    % Convert to constrained zonotope object
    cZ = conZonotope(P);

    % Calculate vertices
    V1 = vertices(cZ);

    % Convert back to a polytope
    P = polytope(cZ);

    % Calculate vertices
    V2 = vertices(P);

    % plot the result
%     plot(cZ,[1,2],'FaceColor','b');
%     hold on
%     plot(P,[1,2],'r');
%     plot(V(1,:),V(2,:),'.k','MarkerSize',12);

    % check for correctness
    % note that tolerance is quite large since computation of vertices uses
    % support functions, i.e., linear programs, and 1e-4 relates to that
    % solver tolerance
    assertLoop(all(contains(P,V1,'exact',1e-4)),j)
    assertLoop(all(contains(cZ,V2,'exact',1e-4)),j)
end

% ------------------------------ END OF CODE ------------------------------
