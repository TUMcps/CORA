function simRes = plus(simRes,v)
% plus - Overloaded '+' operator for the states of the simulated trajectories
%
% Syntax:
%    simRes = plus(simRes,v)
%
% Inputs:
%    simRes - simResult object
%    S - contSet object or numerical vector
%
% Outputs:
%    simRes - transformed simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get reachSet object
[simRes,v] = findClassArg(simRes,v,'simResult');

if ~isnumeric(v)
    throw(CORAerror("CORA:noops",simRes,v))
end

% compute Minkowski sum
for r=1:length(simRes)
    for i=1:length(simRes(r).x)
        simRes(r).x{i} = simRes(r).x{i} + v';
    end
end

end

% ------------------------------ END OF CODE ------------------------------
