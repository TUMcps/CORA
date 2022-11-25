function p0=initialProbability(obj)
% initialProbability - Calculate probability vector out of the 
% state space partition and the initial state set
%
% Syntax:  
%    p0=initialProbability(obj)
%
% Inputs:
%    obj - simulation object
%
% Outputs:
%    p0 - initial probability vector
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 09-October-2006
% Last update: 17-June-2008 
% Last revision: ---

%------------- BEGIN CODE --------------

[~, p0] = exactIntersectingCells(obj.simOptions.stateField, obj.simOptions.initialStateSet);
    
%------------- END OF CODE --------------