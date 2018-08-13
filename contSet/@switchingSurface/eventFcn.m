function [value,isterminal,direction] = eventFcn(obj,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
% trajectory enters or leaves a halfspace;
% this event function is needed, e.g. for matlab ode-solvers
%
% Syntax:  
%    [value,isterminal,direction] = eventFcn(obj,x,direction)
%
% Inputs:
%    obj - halfspace object
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

% Author:       Matthias Althoff
% Written:      06-June-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute value to determine zero crossing
value = obj.c'*x-obj.d;
% Always stop the integration when event detected
isterminal = 1;   
% direction is fed through and not changed

%------------- END OF CODE --------------