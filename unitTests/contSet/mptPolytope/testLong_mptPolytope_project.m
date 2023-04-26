function res = testLong_mptPolytope_project
% testLong_mptPolytope_project - unit test function for projection
%    of polytopes
%
% Syntax:  
%    res = testLong_mptPolytope_project()
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

% Author:       Niklas Kochdumper
% Written:      21-December-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% compare projection with the projection of the polytope vertices and check
% if the same result is obtained
res = true;

for i = 1:2
   
    % create random polytope
    poly = mptPolytope.generateRandom('Dimension',3);
    
    % create random projection dimensions
    dims = randperm(3);
    dims = dims(1:i);
    
    % project using mptPolytope/project function
    poly_ = project(poly,dims);
    
    % project polytope using the polytope vertices
    V = vertices(poly);
    V = V(dims,:);
    
    if i == 1
        V = [min(V),max(V)];
    else
        k = convhulln(V');
        V = V(:,unique(k));
    end
    
    % compute vertices of projection
    V_ = vertices(poly_);
    
    % compare the vertices
    if ~compareMatrices(V,V_,1e-12)
        res = false; break
    end
end

%------------- END OF CODE --------------
