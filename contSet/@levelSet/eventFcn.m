function [value,isterminal,direction] = eventFcn(ls,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
%    trajectory enters or leaves a interval; this event function is needed,
%    e.g., for MATLAB ODE solvers
%
% Syntax:
%    [value,isterminal,direction] = eventFcn(ls,x,direction)
%
% Inputs:
%    ls - levelSet object
%    x - system state
%    direction - event if the state enters or leaves the set
%
% Outputs:
%    value - value of the event function
%    isterminal - specifies if the simulation stops if an event turns zero
%    direction - specifies if the value of the event function has to 
%                turn from negative to positive or the other way round
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       20-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% detect zero crossings
value = ls.funHan(x);

% always stop the integration when event detected
isterminal = ones(length(value),1);   

% vectorize direction
direction = ones(length(value),1)*direction; 

% ------------------------------ END OF CODE ------------------------------
