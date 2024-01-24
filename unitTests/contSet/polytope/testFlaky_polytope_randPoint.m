function res = testFlaky_polytope_randPoint
% testFlaky_polytope_randPoint - unit test function of random point
%    sampling
%
% Syntax:
%    res = testFlaky_polytope_randPoint
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
% Written:       30-November-2022
% Last update:   22-May-2023 (AK, new methods)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% tolerance
tol = 1e-8;

% number of tests
nrTests = 50;

for i=1:nrTests

    % random dimension
    n = randi([2,5]);

    % init polytope
    P = polytope.generateRandom('Dimension',n);

    % sample random points using different syntaxes
    p = [];
    p = randPoint(P);
    p(:,end+1) = randPoint(P,1);
    p(:,end+1) = randPoint(P,1,'standard');
    p(:,end+1:end+3) = randPoint(P,3);
    p(:,end+1) = randPoint(P,1,'extreme');
    p(:,end+1) = randPoint(P,1,'uniform');
    p(:,end+1) = randPoint(P,1,'uniform:hitAndRun');
    p(:,end+1) = randPoint(P,1,'uniform:billiardWalk');

    % all need to be contained within the polytope
    if ~all(contains(P,p,'exact',tol))
        res = false; return
    end
end


% number of tests
nrOfTests = 50;
methods = {'standard','extreme','uniform',...
    'uniform:hitAndRun','uniform:billiardWalk'};

for i=1:length(methods)
    for j=1:nrOfTests
        % random dimension
        n = randi([2,6]); % small because of containment checks
        
        % instantiate random polytope
        P = polytope.generateRandom('Dimension', n);
        
        % compute random points
        nrPts = 10;
        p = randPoint(P,nrPts,methods{i});
        
        % check for containment in polytope
        res = all(contains(P,p,'exact',tol));
        
    end
end

% ------------------------------ END OF CODE ------------------------------
