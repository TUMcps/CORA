function [val,isterminal,direction] = eventFcn(fs,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
%    trajectory enters or leaves a fullspace (can never leave);
%    this event function is needed, e.g. for matlab ODE-solvers
%
% Syntax:
%    [val,isterminal,direction] = eventFcn(fs,x,direction)
%
% Inputs:
%    fs - fullspace object
%    x - system state
%    direction - event if the state enters or leaves the set
%
% Outputs:
%    val - value of the event function
%    isterminal - specifies if the simulation stops if an event turns zero
%    direction - specifies if the value of the event function has to 
%                turn from negative to positive or the other way round
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/eventFcn

% Authors:       Mark Wetzlinger
% Written:       14-December-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute value
val = -Inf(fs.dimension,1);
% always stop the integration when event detected
isterminal = ones(length(val),1);   
% vectorize direction
direction = ones(length(val),1)*direction; 

% ------------------------------ END OF CODE ------------------------------
