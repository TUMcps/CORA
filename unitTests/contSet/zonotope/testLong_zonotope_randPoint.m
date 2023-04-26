function res = testLong_zonotope_randPoint
% testLong_zonotope_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLong_zonotope_randPoint
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

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  14-March-2021 (MW, adapt to new syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% tolerance
tol = 1e-9;

% check empty zonotope object
res = isempty(randPoint(zonotope()));


% number of tests
nrOfTests = 100;

for i=1:nrOfTests
    % random dimension
    n = randi([2,8]); % small because of containment checks
    
    % random center
    c = randn(n,1);
    % random generator matrix -> parallelotope
    G = randn(n);
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % check 'standard' method
    nrPtsStandard = 25;
    pNormal = randPoint(Z,nrPtsStandard,'standard');
    
    % check 'extreme' method:
    % less points than extreme points
    nrPtsExtreme(1) = ceil(2^n * 0.5);
    pExtreme = randPoint(Z,nrPtsExtreme(1),'extreme');
    % as many points as extreme points
    nrPtsExtreme(2) = ceil(2^n);
    pExtreme = [pExtreme, randPoint(Z,nrPtsExtreme(2),'extreme')];
    % more points than extreme points
    nrPtsExtreme(3) = ceil(2^n * 5);
    pExtreme = [pExtreme, randPoint(Z,nrPtsExtreme(3),'extreme')];
    nrPtsExtreme = sum(nrPtsExtreme);
    
    % check 'gaussian' method:
    nrPtsGaussian = 25;
    pr = 0.8;
    pGaussian = randPoint(Z,nrPtsGaussian,'gaussian',pr);
    % note: these points are not checked for containment because there
    %       are not guaranteed to be inside of Z
    
    % check for containment in zonotope
    % (use containsPoint instead of in, since the latter has a bug)
    if ~all(contains(Z,pNormal))
        res = false; break;
    end
    if ~all(contains(enlarge(Z,1+tol),pExtreme))
        % enlarging Z is required, otherwise wrong result!
        res = false; break;
    end

end

%------------- END OF CODE --------------
