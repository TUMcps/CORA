function res = testLong_polyZonotope_randPoint
% testLong_polyZonotope_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testLong_polyZonotope_randPoint
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

% Authors:       Mark Wetzlinger
% Written:       02-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-9;

% check empty zonotope object -> error
res = isempty(randPoint(polyZonotope.empty(1)));


% number of tests
nrOfTests = 100;

for i=1:nrOfTests
    % random dimension
    n = randi([2,5]);
    % random number of generators
    nrIndepGens = randi(25);
    nrDepGens = randi(25);
    % random number of dependent factors
    nrDepFactors = randi([nrDepGens,nrDepGens+10]);
    
    % random center, random generator matrices
    c = randn(n,1);
    G = randn(n,nrDepGens);
    GI = randn(n,nrIndepGens);
    % random exponent matrix (no redudancies) and indentifiers
    counter = 0;
    while true
        counter = counter + 1;
        E = randi(7,nrDepFactors,nrDepGens) - 2;
        E(E < 0) = 0;
        if size(unique(E','rows')',2) == size(E,2)
            break;
        end
        if counter > 10
            continue; % skip current n/nrDepGens/nrDepFactors combination
        end
    end
    id = randperm(nrDepFactors)';
    
    % center, dependent, independent generators, and exponent matrix
    pZ = polyZonotope(c,G,GI,E);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(pZ.GI,GI,tol) ...
            || ~compareMatrices([G;E],[pZ.G;pZ.E],tol)
        res = false; return;
    end
    
    
    % check 'standard' method
    nrPtsStandard = 100;
    pNormal = randPoint(pZ,nrPtsStandard,'standard');
    
    % check 'extreme' method:
    % less points than extreme points
    nrPtsExtreme(1) = ceil(2^n * 0.5);
    pExtreme = randPoint(pZ,nrPtsExtreme(1),'extreme');
    % as many points as extreme points
    nrPtsExtreme(2) = ceil(2^n);
    pExtreme = [pExtreme, randPoint(pZ,nrPtsExtreme(2),'extreme')];
    % more points than extreme points
    nrPtsExtreme(3) = ceil(2^n * 5);
    pExtreme = [pExtreme, randPoint(pZ,nrPtsExtreme(3),'extreme')];
    nrPtsExtreme = sum(nrPtsExtreme);

    
    % check for containment in polyZonotope
    % NOTE: add this to the test once point in pZ is implemented
%     for j=1:nrPtsStandard
%         if ~contains(pZ,pNormal(:,j))
%             res = false; return
%         end
%     end
%     for j=1:nrPtsExtreme
%         if ~contains(enlarge(pZ,1+tol),pExtreme(:,j))
%             % enlarging Z is required, otherwise wrong result!
%             res = false; return
%         end
%     end
end

% ------------------------------ END OF CODE ------------------------------
