function [inputProb]=priv_inputDist(p)
% priv_inputDist - Returns the input probability distribution of the simulated
% traffic participant
%
% Syntax:
%    [inputProb]=priv_inputDist(p)
%
% Inputs:
%    p - joint probabilities (state and input) of different time steps
%
% Outputs:
%    inputProb - probability distribution of the input
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       24-April-2009
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for iStep=1:length(p)
    %project
    inputProb{iStep}=sum(p{iStep});
end

% ------------------------------ END OF CODE ------------------------------
