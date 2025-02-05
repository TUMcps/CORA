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
% See also: none

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       17-September-2019
% Last update:   14-March-2021 (MW, adapt to new syntax)
%                22-May-2023 (AK, added uniform sampling)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-6;

% number of tests
nrOfTests = 10;

for i=1:nrOfTests

    % random dimension
    n = randi([2,8]); % small because of containment checks

    % random zonotope
    Z = zonotope.generateRandom("Dimension",n,"NrGenerators",2*n);
    
    % check 'standard' method
    nrPtsStandard = 25;
    pStandard = randPoint(Z,nrPtsStandard,'standard');
    
    % check 'uniform' methods
    nrPtsUniform = 10;
    pUniform = randPoint(Z,nrPtsUniform,'uniform');
    pUniformHitAndRun = randPoint(Z,nrPtsUniform,'uniform:hitAndRun');
    pUniformBilliard = randPoint(Z,nrPtsUniform,'uniform:billiardWalk');
    pUniformBallWalk = randPoint(Z,nrPtsUniform,'uniform:ballWalk');
    
    
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
    % note: these points are not checked for containment because they
    %       are not guaranteed to be inside of Z
    
    
    % check for containment in zonotope
    assertLoop(all(contains(Z,pStandard)),i)
    assertLoop(all(contains(Z,pUniform)),i)
    assertLoop(all(contains(Z,pUniformHitAndRun)),i)
    assertLoop(all(contains(Z,pUniformBilliard)),i)
    assertLoop(all(contains(Z,pUniformBallWalk)),i)
    assertLoop(all(contains(Z,pExtreme,'exact',tol)),i)

end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
