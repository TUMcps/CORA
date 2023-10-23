function res = testLong_polytope_project
% testLong_polytope_project - unit test function for projection
%    of polytopes
%
% Syntax:
%    res = testLong_polytope_project()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       21-December-2020
% Last update:   05-December-2022 (MW, integrate comparison to zonotope)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

% idea: compare vertices of projection with projected vertices
for i=1:nrTests
    
    % random dimension (low... to facilitate vertex enumeration)
    n = randi([3,4]);

    % create random polytope
    P = polytope.generateRandom('Dimension',n);
    
    % create random projection dimensions
    projDims = randperm(n);
    % either one or two dimensions removed
    projDims = projDims(1:randi(2));
    
    % computation of vertices or projection
    % project using polytope/project function
    P_ = project(P,projDims);
    % compute vertices of projection
    V_ = vertices(P_);
    
    % computation of projected vertices
    V = vertices(P);
    V = V(projDims,:);
    % convex hull to remove redundancies
    if length(projDims) == 1
        V = [min(V),max(V)];
    else
        k = convhulln(V');
        V = V(:,unique(k));
    end
    
    % compare vertices
    if ~compareMatrices(V_,V,1e-12)
        res = false; break
    end


    % compare to interval projection
    I = interval.generateRandom('Dimension',n);
    % convert to polytope
    P = polytope(I);

    % project both sets
    I_ = project(I,projDims);
    P_ = project(P,projDims);

    % compute vertices
    V_I = vertices(I_);
    V_P = vertices(P_);

    % compare vertices
    if ~compareMatrices(V_I,V_P,1e-12)
        res = false; break
    end


    % compare to zonotope projection
    Z = interval.generateRandom('Dimension',n);
    % convert to polytope
    P = polytope(Z);

    % project both sets
    Z_ = project(Z,projDims);
    P_ = project(P,projDims);

    % compute vertices
    V_Z = vertices(Z_);
    V_P = vertices(P_);

    % compare vertices
    if ~compareMatrices(V_Z,V_P,1e-12)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
