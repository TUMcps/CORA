function res = testLong_zonotope_project
% testLong_zonotope_project - unit test function of project
%
% Syntax:
%    res = testLong_zonotope_project
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

for i=1:nrTests

    % random dimension
    n = randi([2,10]);

    % create a random zonotope
    nrOfGens = randi([2*n,5*n]);
    c = -1+2*rand(n,1);
    G = -1+2*rand(n,nrOfGens);
    Z = zonotope(c,G);

    % choose random subspace
    projDims = randperm(n,min([n-1,2]));
    if rand < 0.5
        % choose logical indexing
        temp = projDims;
        projDims = false(n,1);
        projDims(temp) = true;
    end

    % project original center and generator matrix
    cproj = c(projDims);
    Gproj = G(projDims,:);
    
    % project zonotope
    Zproj = project(Z,projDims);
    
    % check projections return the same result
    if ~compareMatrices(cproj,center(Zproj)) ...
            || ~compareMatrices(Gproj,generators(Zproj))
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
