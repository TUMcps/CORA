function res = testLong_conZonotope_mptPolytope
% testLong_conZonotope_mptPolytope - unit test function for
%    conversion between constrained zonotopes and polytopes
%
% Syntax:  
%    res = testLong_conZonotope_mptPolytope
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

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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

    % Construct a mptPolytope object from the vertices
    P = mptPolytope(V');

    % Convert to constrained zonotope object
    cZ = conZonotope(P);

    % Calculate vertices
    V1 = vertices(cZ);

    % Convert back to a mptPolytope
    P = mptPolytope(cZ);

    % Calculate vertices
    V2 = vertices(P);

    % plot the result
%     plot(cZ,[1,2],'FaceColor','b');
%     hold on
%     plot(P,[1,2],'r');
%     plot(V(1,:),V(2,:),'.k','MarkerSize',12);

    % check for correctness
    % note that tolerance is quite large since computation of vertices uses
    % support functions, i.e., linear programs, and 1e-6 relates to that
    % solver tolerance
    if n == 3
        if ~compareMatrices(V,V1,1e-6)
            throw(CORAerror('CORA:testFailed'));
        elseif ~compareMatrices(V,V2,1e-6)
            throw(CORAerror('CORA:testFailed'));
        end
    elseif n == 2
        % 2D check for correctness using polygons (more robust handling of
        % collinear points, see projVertices)
        poly = polygon(V); poly1 = polygon(V1); poly2 = polygon(V2);
        if ~isequal(poly,poly1,1e-6)
            throw(CORAerror('CORA:testFailed'));
        elseif ~isequal(poly,poly2,1e-6)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

%------------- END OF CODE --------------