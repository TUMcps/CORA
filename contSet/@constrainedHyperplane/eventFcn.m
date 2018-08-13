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
% Written:      10-August-2011
% Last update:  14-May-2017
% Last revision:---

%------------- BEGIN CODE --------------

%extract hyperplane values
HP_c = obj.h.c;
HP_d = obj.h.d;

%compute value to determine zero crossing
value = HP_c'*x - HP_d;

isterminal = 1;

% Old setting (changed due to changed simulation technique):
% %compute constraint satisfaction
% if all(obj.C*x - obj.d <= 0)
%     isterminal = 1;
% else
%     isterminal = 0;
% end


%------------- END OF CODE --------------