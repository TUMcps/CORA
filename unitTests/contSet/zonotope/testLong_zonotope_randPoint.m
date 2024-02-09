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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       17-September-2019
% Last update:   14-March-2021 (MW, adapt to new syntax)
%                22-May-2023 (AK, added uniform sampling)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-6;

% check empty zonotope object
res = isempty(randPoint(zonotope.empty(1)));


% number of tests
nrOfTests = 100;

for i=1:nrOfTests

    % random dimension
    n = randi([2,8]); % small because of containment checks

    % random center
    c = randn(n,1);
    % random generator matrix -> parallelotope
    G = randn(n);
    while abs(det(G)) < 1e-4
        G = randn(n);
    end
    
    % instantiate zonotope
    % For the first case, we try with an empty zonotope; for the second, we
    % try with one that has no generators. For the third, we try with a
    % zonotope that is degenerate. The fourth is a parallelotope, the rest
    % are zonotopes that are not parallelotopes.
    if i == 1
        Z = zonotope.empty(1);
    elseif i == 2
        Z = zonotope(c, []);
    elseif i == 3
        Z = zonotope([c;0], [G;zeros([1 size(G,2)])]);
    elseif i == 4
        Z = zonotope(c,G);
    else
        Grest = randn(n);
        Z = zonotope(c,[G Grest]);
    end
    
    % check 'standard' method
    nrPtsStandard = 25;
    pStandard = randPoint(Z,nrPtsStandard,'standard');
    
    % check 'uniform' methods
    nrPtsUniform = 10;
    pUniform = randPoint(Z,nrPtsUniform,'uniform');
    pUniformHitAndRun = randPoint(Z,nrPtsUniform,'uniform:hitAndRun');
    pUniformBilliard = randPoint(Z,nrPtsUniform,'uniform:billiardWalk');
    
    
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
    if i ~= 3
        % TODO: Fix the Gaussian sampling, so that it also works for
        % degenerate sets
        pGaussian = randPoint(Z,nrPtsGaussian,'gaussian',pr);
    end
    % note: these points are not checked for containment because they
    %       are not guaranteed to be inside of Z
    
    
    % check for containment in zonotope
    if ~all(contains(Z,pStandard))
        res = false; break;
    end
    if ~all(contains(Z,pUniform))
        res = false;break;
    end
    if ~all(contains(Z,pUniformHitAndRun))
        res = false; break;
    end
    if ~all(contains(Z,pUniformBilliard))
        res = false; break;
    end
    if ~all(contains(Z,pExtreme,'exact',tol))
        % enlarging Z is required, otherwise wrong result!
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
