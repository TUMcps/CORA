function simRes = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the states of the simulated trajectories
%
% Syntax:
%    simRes = mtimes(A,simRes)
%
% Inputs:
%    A - numeric matrix
%    simRes - simResult object
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

[simRes,A] = findClassArg(factor1,factor2,'simResult');

if ~isnumeric(A)
    throw(CORAerror("CORA:noops",A,simRes))
end

% compute linear map
for r=1:length(simRes)
    for i=1:length(simRes(r).x)
        simRes(r).x{i} = simRes(r).x{i} * A';
    end
end

end

% ------------------------------ END OF CODE ------------------------------
