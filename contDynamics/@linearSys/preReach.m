function [obj] = preReach(obj,options)
% preReach - prepares reachable set computation for linear systems
%
% Syntax:  
%    [obj] = preReach(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-May-2011
% Last update:  29-October-2017
% Last revision:---

%------------- BEGIN CODE --------------


% compute exponential matrix
obj = exponential(obj,options);
% compute time interval error (tie)
obj = tie(obj,options);
% compute reachable set due to input
obj = inputSolution(obj,options);
%change the time step
obj.taylor.timeStep=options.timeStep;
% compute exponential matrix
obj.taylor.eAt=expm(obj.A*options.timeStep);


%------------- END OF CODE --------------