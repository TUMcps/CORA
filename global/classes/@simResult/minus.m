function simRes = minus(simRes1,simRes2)
% minus - Overloaded '-' operator for the states of the simulated trajectories
%
% Syntax:
%    simRes = minus(simRes,factor2)
%
% Inputs:
%    simRes1 - numeric or simResult object
%    factor2 - numeric or simResult object 
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

simRes = simRes1 + (-simRes2);

end

% ------------------------------ END OF CODE ------------------------------
