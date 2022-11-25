function [value,isterminal,direction] = eventFcn(obj,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
% trajectory enters or leaves a interval;
% this event function is needed, e.g. for matlab ode-solvers
%
% Syntax:  
%    [value,isterminal,direction] = eventFcn(obj,x,direction)
%
% Inputs:
%    obj - levelSet object
%    x - system state
%    direction - event if the state enters or leaves the set
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

% Author:       Niklas Kochdumper
% Written:      20-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % detect zero crossings
    value = obj.funHan(x);

    % always stop the integration when event detected
    isterminal = ones(length(value),1);   

    % vectorize direction
    direction = ones(length(value),1)*direction; 

%------------- END OF CODE --------------