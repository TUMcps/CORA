function res = testLongDuration_polyZonotope_randPoint
% testLongDuration_polyZonotope_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLongDuration_polyZonotope_randPoint
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
% Written:      02-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% tolerance
tol = 1e-9;

% check empty zonotope object -> error
res = isempty(randPoint(polyZonotope()));


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
    Grest = randn(n,nrIndepGens);
    % random exponent matrix (no redudancies) and indentifiers
    counter = 0;
    while true
        counter = counter + 1;
        expMat = randi(7,nrDepFactors,nrDepGens) - 2;
        expMat(expMat < 0) = 0;
        if size(unique(expMat','rows')',2) == size(expMat,2)
            break;
        end
        if counter > 10
            continue; % skip current n/nrDepGens/nrDepFactors combination
        end
    end
    id = randperm(nrDepFactors)';
    
    % center, dependent, independent generators, and exponent matrix
    pZ = polyZonotope(c,G,Grest,expMat);
    if ~all(withinTol(center(pZ),c,tol)) || ~compareMatrices(pZ.Grest,Grest,tol) ...
            || ~compareMatrices([G;expMat],[pZ.G;pZ.expMat],tol)
        res = false; break;
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
%             res = false; break;
%         end
%     end
%     for j=1:nrPtsExtreme
%         if ~contains(enlarge(pZ,1+tol),pExtreme(:,j))
%             % enlarging Z is required, otherwise wrong result!
%             res = false; break;
%         end
%     end
%     
%     if ~res
%         break;
%     end
end

if ~res
    path = pathFailedTests(mfilename());
    save(path,'pZ','Grest','c');
end

%------------- END OF CODE --------------
