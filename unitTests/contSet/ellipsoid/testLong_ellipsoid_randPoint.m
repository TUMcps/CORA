function res = testLong_ellipsoid_randPoint
% testLong_ellipsoid_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testLong_ellipsoid_randPoint
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

% Authors:       Victor Gassmann
% Written:       02-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true;

% number of tests
nrOfTests = 100;
nrPts = 10;

for i=1:nrOfTests
    %%% generate all variables necessary to replicate results
    % random dimension
    n = randi([2,8]); % small because of containment checks
    % random ellipsoid
    E = ellipsoid.generateRandom('Dimension',n);
    % random ellipsoid that is degenerate
    Ed = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',true); 
    % compute random points
    pNormal = randPoint(E,nrPts,'standard');
    pExtreme = randPoint(E,nrPts,'extreme');
    pNormal_d = randPoint(Ed,nrPts,'standard');
    pExtreme_d = randPoint(Ed,nrPts,'extreme');
    
    % check for containment in ellipsoid
    if ~all(contains(E,pNormal)) || ~all(contains(Ed,pNormal_d)) || ...
       ~all(contains(enlarge(E,1+tol),pExtreme)) || ~all(contains(enlarge(Ed,1+tol),pExtreme_d))
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
