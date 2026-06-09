function res = test_zonoBundle_projVertices
% test_zonoBundle_projVertices - unit test function for computation of
%    vertices of a 2D projection
%
% Syntax:
%    res = test_zonoBundle_projVertices
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
    
    % 2D zonotope bundle
    Z1 = zonotope([1 3 0 1; 1 0 2 1]);
    Z2 = zonotope([0 2 2; 0 2 -2]);
    
    zB = zonoBundle({Z1,Z2});
    
    % compute vertices
    V = vertices(zB);
    V_proj1 = projVertices(zB,[1,2],'angle');
    V_proj2 = projVertices(zB,[1,2],'supportFunc');
    
    % check vertices
    assert(compareMatrices(V,V_proj2,1e-14,'subset'))
    assert(compareMatrices(V,V_proj1,1e-14,'subset'))
    
    
    % higher-dimensional zonotope bundle
    Z1 = zonotope([1 3 0 1 0; 1 0 2 1 0; 1 0 0 0 1]);
    Z2 = zonotope([0 2 2 0; 0 2 -2 0; 2 0 0 0.5]);
    
    zB = zonoBundle({Z1,Z2});
    
    % compute vertices of full zonotope
    V = vertices(zB);
    
    % dimensions for projection
    dims = {[1,2],[2,3],[1,3]};
    
    % check all three projections
    for i = 1:length(dims)
    
        % computed vertices of projected constrained zonotope
        V_proj1 = projVertices(zB,dims{i},'angle');
        V_proj2 = projVertices(zB,dims{i},'supportFunc');
    
        % check vertices
        assertLoop(all(contains(polytope(V_proj1),V(dims{i},:), ...
                                                        'exact',1e-6)),i);
        assertLoop(all(contains(polytope(V(dims{i},:)),V_proj1, ...
                                                        'exact',1e-6)),i);
        assertLoop(all(contains(polytope(V_proj2),V(dims{i},:), ...
                                                        'exact',1e-6)),i);
        assertLoop(all(contains(polytope(V(dims{i},:)),V_proj2, ...
                                                        'exact',1e-6)),i);
    end

end

% ------------------------------ END OF CODE ------------------------------
