function res = testFlaky_polytope_project
% testFlaky_polytope_project - unit test function for projection
%    of polytopes
%
% Syntax:
%    res = testFlaky_polytope_project()
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
%                21-November-2023 (MW, replace vertex computation by support function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

% number of support function evaluations
nrsF = 10;

% idea: compare vertices of projection with projected vertices
for i=1:nrTests
    
    % low random dimension to facilitate vertex enumeration
    n = randi([3,4]);

    % create random polytope
    P = polytope.generateRandom('Dimension',n);
    
    % create random projection dimensions
    projDims = randperm(n);
    % either one or two dimensions removed
    projDims = projDims(1:randi(2));
    
    % project using polytope/project function
    P_ = project(P,projDims);

    % evaluate support function of original set and projected set in
    % directions that lie in the projected subspace
    for j=1:nrsF
        % random normalized direction in subspace
        l = randn(n,1);
        l(~any((1:n)'==projDims,2)) = 0;
        l = l./vecnorm(l);
        l_ = l(projDims);

        % compare values support function
        val = supportFunc(P,l);
        val_ = supportFunc(P_,l_);
        if ~withinTol(val,val_,1e-5)
            throw(CORAerror('CORA:testFailed'));
        end
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
        throw(CORAerror('CORA:testFailed'));
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
        throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
