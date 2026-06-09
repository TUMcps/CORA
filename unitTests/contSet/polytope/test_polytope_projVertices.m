function res = test_polytope_projVertices
% test_polytope_projVertices - unit test function for computation of
%    vertices of a 2D projection
%
% Syntax:
%    res = test_polytope_projVertices
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
% See also: none

% Authors:       Niklas Kochdumper
% Written:       20-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % assume true
    res = true;


    % polytope with equality constraint
    % higher-dimensional constrained zonotope
    A = [1 1 1;-1 0 0;0 -1 0;0 0 -1];
    b = [1;0;0;0];

    Aeq = [0.1 -0.1 1];
    beq = 0.5;
    
    P = polytope(A,b,Aeq,beq);

    % compute vertices
    V = vertices(P);
    V_proj = projVertices(P);
    
    % check vertices
    assert(compareMatrices(V_proj,V(1:2,:),1e-14,'subset'))

    
    % higher-dimensional constrained zonotope
    A = [1 1 1;-1 0 0;0 -1 0;0 0 -1];
    b = [1;0;0;0];
    
    P = polytope(A,b);
    
    % compute vertices of full constrained zonotope
    V = vertices(P);
    
    % dimensions for projection
    dims = {[1,2],[2,3],[1,3]};
    
    % check all three projections
    for i = 1:length(dims)
    
        % computed vertices of projected polytope
        V_proj = projVertices(P,dims{i});
    
        % check vertices
        assertLoop(compareMatrices(V_proj,V(dims{i},:),1e-6,'subset'),i)
    end
    
    
    % convex hull of two polytopes
    A = [1 1 1;-1 0 0;0 -1 0;0 0 -1];
    b = [1;0;0;0];
    P1 = polytope(A,b);
    
    P2 = P1 + [1;2;-1];
    
    % compute convex hull
    P_ = convHull(P1,P2);
    
    % compute vertices
    V = vertices(P_);
    V_proj = projVertices(P_);
    
    % check vertices
    assert(compareMatrices(V_proj,V(1:2,:),1e-14,'subset'))

end

% ------------------------------ END OF CODE ------------------------------
