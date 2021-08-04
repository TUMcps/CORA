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
%    res - boolean 
%
% Example: 
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

res = true;
tol = 1e-9;

% check empty zonotope object -> error
pZ = polyZonotope();
try 
    p = randPoint(pZ); % <- should throw error here
    res = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        res = false;
    end
end


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
    if any(abs(center(pZ) - c) > tol) || any(any(abs(pZ.Grest - Grest) > tol)) ...
            || ~sameGexpMat(G,expMat,pZ.G,pZ.expMat,tol)
        res_rand = false; break;
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
%         if ~in(pZ,pNormal(:,j))
%             res = false; break;
%         end
%     end
%     for j=1:nrPtsExtreme
%         if ~in(enlarge(pZ,1+tol),pExtreme(:,j))
%             % enlarging Z is required, otherwise wrong result!
%             res = false; break;
%         end
%     end
%     
%     if ~res
%         break;
%     end
end


if res
    disp('testLongDuration_polyZonotope_randPoint successful');
else
    disp('testLongDuration_polyZonotope_randPoint failed');
end

end


% Auxiliary Function ------------------------------------------------------
function same = sameGexpMat(G,expMat,G_new,expMat_new,tol)
% assumption: no redundancies removal in expMat -> expMat_new

if any(size(G) ~= size(G_new)) || any(size(expMat) ~= size(expMat_new))
    same = false; return;
end

% quick check (same order as before)
if ~(any(any(abs(G_new - G) > tol)) || any(any(abs(expMat_new - expMat) > tol)))
    same = true; return;
end

% init output argument
same = true;
% go through generators and exponent columns one by one
nrGens = size(G,2);
for k=1:nrGens
    % select first generator
    gk = G(:,1);
    
    % search for gk in G_new
    idx = find(all(abs(G_new - gk) < tol,1),1,'first');
    % generator not found
    if isempty(idx)
        same = false; break;
    end
    
    % check equality of expMat columns
    if any(abs(expMat(:,1) - expMat_new(:,idx)) > tol)
        same = false; break;
    end
    
    % remove columns in Gs and expMats
    G(:,1) = []; expMat(:,1) = [];
    G_new(:,idx) = []; expMat_new(:,idx) = [];
end

end

%------------- END OF CODE --------------
