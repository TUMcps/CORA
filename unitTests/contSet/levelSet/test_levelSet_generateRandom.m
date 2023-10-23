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

resvec = [];

try 
    % empty call
    ls = levelSet.generateRandom();
    resvec(end+1)=true;
    
    % check dimension
    n = 3;
    ls = levelSet.generateRandom('Dimension',n);
    resvec(end+1) = dim(ls) == n;
    
    % check nrEquations
    nrEqs = 4;
    ls = levelSet.generateRandom('NrEquations',nrEqs);
    resvec(end+1) = length(ls.eq) == nrEqs;
    
    % check CompOps
    compOps = {'<'};
    ls = levelSet.generateRandom('CompOps',compOps);
    resvec(end+1) = all(cellfun(@(x) ismember(x,compOps),ls.compOp,'UniformOutput',true));

    % check combination
    ls = levelSet.generateRandom('Dimension',n,'NrEquations',nrEqs,'CompOps',compOps);
    resvec(end+1) = dim(ls) == n;
    resvec(end+1) = length(ls.eq) == nrEqs;
    resvec(end+1) = all(cellfun(@(x) ismember(x,compOps),ls.compOp,'UniformOutput',true));

catch ME
    resvec(end+1) = false;
end

% unify results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
