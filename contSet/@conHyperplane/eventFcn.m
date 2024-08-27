function [val,isterminal,direction] = eventFcn(hyp,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
%    trajectory crosses a constrained hyperplane;
%    this event function is needed, e.g. for matlab ODE-solvers
%
% Syntax:
%    [val,isterminal,direction] = eventFcn(P,x,direction)
%
% Inputs:
%    hyp - conHyperplane object
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

% Authors:       Niklas Kochdumper
% Written:       25-July-2024 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % compute distance to hyperplane
    val = hyp.a*x - hyp.b;

    % always stop the integration when event detected
    isterminal = 1;   

% ------------------------------ END OF CODE ------------------------------
