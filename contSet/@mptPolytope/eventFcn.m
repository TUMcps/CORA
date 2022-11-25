function [value,isterminal,direction] = eventFcn(obj,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
% trajectory enters or leaves a mptPolytope;
% this event function is needed, e.g. for matlab ode-solvers
%
% Syntax:  
%    [value,isterminal,direction] = eventFcn(obj,x,direction)
%
% Inputs:
%    obj - mptPolytope object
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

% Author:       Matthias Althoff
% Written:      12-February-2012 
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

try %MPT V3
    H = obj.P.A;
    K = obj.P.b;
catch %MPT V2
    [H,K] = double(obj.P);
end

value=H*x-K;
% Always stop the integration when event detected
isterminal = ones(length(K),1);   
% Vectorize direction
direction = ones(length(K),1)*direction; 

%------------- END OF CODE --------------
