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

    % Check for correctness
    if ~compareMatrices(V,V1,1e-10)
       throw(CORAerror('CORA:testFailed'));
    elseif ~compareMatrices(V,V2,1e-10)
       throw(CORAerror('CORA:testFailed'));
    end
    
end

%------------- END OF CODE --------------