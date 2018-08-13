function [handle] = eventFcn(obj,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
% trajectory enters or leaves a zonotope;
% this event function is needed, e.g. for matlab ode-solvers
%
% Syntax:  
%    [value,isterminal,direction] = eventFcn(obj,x,direction)
%
% Inputs:
%    obj - zonotope object
%    x - system state
%    direction - event if the state enters or leaves the interval hull?
%
% Outputs:
%    value - value of the event function
%    isterminal - specifies if the simulation stops if an event turns zero
%    direction - specifies if the value of the event function has to 
%    turn from negative to positive or the other way round
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

value=obj.halfspace.H*x-obj.halfspace.K;
% Always stop the integration when event detected
isterminal = ones(length(obj.halfspace.K),1)*all(value<=0);   
% Vectorize direction
direction = ones(length(obj.halfspace.K),1)*direction;  

%------------- END OF CODE --------------