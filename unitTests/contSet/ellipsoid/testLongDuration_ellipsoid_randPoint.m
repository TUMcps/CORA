function res = testLongDuration_ellipsoid_randPoint
% testLongDuration_ellipsoid_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLongDuration_ellipsoid_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      02-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;

% check empty ellipsoid object -> empty
res = isempty(randPoint(ellipsoid()));


% number of tests
nrOfTests = 100;

for i=1:nrOfTests
    % random dimension
    n = randi([2,8]); % small because of containment checks
    
    % random ellipsoid
    E = ellipsoid.generateRandom(n);
    
    % random ellipsoid that is degenerate
    Ed = ellipsoid.generateRandom(n,true); 
    
    % compute random points
    nrPts = 10;
    pNormal = randPoint(E,nrPts,'standard');
    pExtreme = randPoint(E,nrPts,'extreme');
    
    pNormal_d = randPoint(Ed,nrPts,'standard');
    pExtreme_d = randPoint(Ed,nrPts,'extreme');
    
    % check for containment in ellipsoid
    if ~in(E,pNormal) || ~in(Ed,pNormal_d) || ...
       ~in(enlarge(E,1+tol),pExtreme) || ~in(enlarge(Ed,1+tol),pExtreme_d)
        res = false;
        break;
    end
end


if res
    disp('test_zonotope_randPoint successful');
else
    disp('test_zonotope_randPoint failed');
end

%------------- END OF CODE --------------
