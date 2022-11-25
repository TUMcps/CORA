function [value,isterminal,direction] = eventFcn(obj,x)
% eventFcn - returns the event function results of a guard set of a 
% transition
%
% Syntax:  
%    [value,isterminal,direction] = eventFcn(obj,x)
%
% Inputs:
%    obj - transition object
%    x - system state
%
% Outputs:
%    value - value of the event function
%    isterminal - specifies if the simulation stops if an event turns zero
%    direction - specifies if the value of the event function has to 
%    turn from negative to positive or the other way round
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: 07-September-2007
% Last revision: ---

%------------- BEGIN CODE --------------

% [value,isterminal,direction] = eventFcn(obj.guard,x,0);
%[value,isterminal,direction] = eventFcn(obj.guard,x,1);
[value,isterminal,direction] = eventFcn(obj.guard,x,-1);

%------------- END OF CODE --------------