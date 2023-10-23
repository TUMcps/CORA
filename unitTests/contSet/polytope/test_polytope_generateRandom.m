function res = test_polytope_generateRandom
% test_polytope_generateRandom - unit test function of random generation;
%    also checks hidden properties and re-instantiates generated polytopes
%    to see if properties are actually correct
%
% Syntax:
%    res = test_polytope_generateRandom
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

% Authors:       Mark Wetzlinger
% Written:       30-November-2022
% Last update:   12-December-2022 (MW, add emptiness check)
%                01-August-2023 (MW, check properties and reinstantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true(0);

% dimension
P = polytope.generateRandom('Dimension',3);
res(end+1,1) = dim(P) == 3 && ~representsa(P,'emptySet');
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet');

% number of constraints
P = polytope.generateRandom('NrConstraints',5);
res(end+1,1) = size(P.A,1) == 5 && ~representsa(P,'emptySet');
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet');

% dimension and number of constraints (unbounded)
P = polytope.generateRandom('Dimension',3,'NrConstraints',3);
res(end+1,1) = dim(P) == 3 && size(P.A,1) == 3 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_);

% dimension and number of constraints (bounded)
P = polytope.generateRandom('Dimension',3,'NrConstraints',10);
res(end+1,1) = dim(P) == 3 && size(P.A,1) == 10 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && P.bounded.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && isBounded(P_);

% dimension and boundedness
P = polytope.generateRandom('Dimension',2,'isBounded',false);
res(end+1,1) = dim(P) == 2 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_);

% dimension, number of constraints, and boundedness
P = polytope.generateRandom('Dimension',3,'NrConstraints',8,'isBounded',false);
res(end+1,1) = dim(P) == 3 && size(P.A,1) == 8 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_);

% degeneracy
P = polytope.generateRandom('IsDegenerate',true);
res(end+1,1) = ~isFullDim(P) ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.fullDim.val) && ~P.fullDim.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isFullDim(P_);

% degeneracy and boundedness
P = polytope.generateRandom('IsBounded',false,'IsDegenerate',true);
res(end+1,1) = ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val ...
    && ~isempty(P.fullDim.val) && ~P.fullDim.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_) && ~isFullDim(P_);

% dimension, degeneracy, and boundedness
P = polytope.generateRandom('Dimension',4,'IsBounded',false,'IsDegenerate',true);
res(end+1,1) = dim(P) == 4 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val ...
    && ~isempty(P.fullDim.val) && ~P.fullDim.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_) && ~isFullDim(P_);

% dimension, number of constraints, degeneracy, and boundedness
P = polytope.generateRandom('Dimension',4,'NrConstraints',10,...
    'IsBounded',false,'IsDegenerate',true);
res(end+1,1) = dim(P) == 4 && size(P.A,1) == 10 ...
    && ~isempty(P.emptySet.val) && ~P.emptySet.val ...
    && ~isempty(P.bounded.val) && ~P.bounded.val ...
    && ~isempty(P.fullDim.val) && ~P.fullDim.val;
P_ = polytope(P.A,P.b,P.Ae,P.be);
res(end+1,1) = ~representsa(P_,'emptySet') && ~isBounded(P_) && ~isFullDim(P_);

% combine results
res = all(res);

% wrong calls, i.e., user inputs cannot be fulfilled

% cannot be degenerate if nrCon = 1
try
    P = polytope.generateRandom('NrConstraints',1,'IsDegenerate',true);
    res = false;
end

% cannot be bounded if nrCon = 1
try
    P = polytope.generateRandom('NrConstraints',1,'IsBounded',true);
    res = false;
end

% cannot be bounded if nrCon < n+1
try
    P = polytope.generateRandom('Dimension',5,'NrConstraints',5,...
        'IsBounded',true);
    res = false;
end

% cannot be bounded and degenerate if nrCon < n+2
try
    P = polytope.generateRandom('Dimension',5,'NrConstraints',6,...
        'IsBounded',true,'IsDegenerate',true);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
