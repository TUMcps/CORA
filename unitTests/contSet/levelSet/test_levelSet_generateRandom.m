function res = test_levelSet_generateRandom
% test_levelSet_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_levelSet_generateRandom
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

% Authors:       Tobias Ladner
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty call
ls = levelSet.generateRandom();

% check dimension
n = 3;
ls = levelSet.generateRandom('Dimension',n);
assert(dim(ls) == n);

% check nrEquations
nrEqs = 4;
ls = levelSet.generateRandom('NrEquations',nrEqs);
assert(length(ls.eq) == nrEqs);

% check CompOps
compOps = {'<'};
ls = levelSet.generateRandom('CompOps',compOps);
assert(all(cellfun(@(x) ismember(x,compOps),ls.compOp,'UniformOutput',true)));

% check combination
ls = levelSet.generateRandom('Dimension',n,'NrEquations',nrEqs,'CompOps',compOps);
assert(dim(ls) == n);
assert(length(ls.eq) == nrEqs);
assert(all(cellfun(@(x) ismember(x,compOps),ls.compOp,'UniformOutput',true)));

% unify results
res = true;

% ------------------------------ END OF CODE ------------------------------
