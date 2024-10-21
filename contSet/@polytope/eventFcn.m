function [val,isterminal,direction] = eventFcn(P,x,direction)
% eventFcn - Returns the results of an event function that detects if a 
%    trajectory enters or leaves a polytope;
%    this event function is needed, e.g. for matlab ODE-solvers
%
% Syntax:
%    [val,isterminal,direction] = eventFcn(P,x,direction)
%
% Inputs:
%    P - polytope object
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
% See also: none

% Authors:       Matthias Althoff
% Written:       12-February-2012 
% Last update:   30-July-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

constraints(P);

if isempty(P.be_.val)
    % compute distance to boundary of each inequality constraint
    val = P.A_.val*x - P.b_.val;
    % Always stop the integration when event detected
    isterminal = ones(length(P.b_.val),1);   
    % Vectorize direction
    direction = ones(length(P.b_.val),1)*direction;

elseif representsa_(P,'conHyperplane',1e-12)
    val = P.Ae_.val*x - P.be_.val;
    isterminal = 1;
    % direction?
    
else
    throw(CORAerror('CORA:notSupported','Given polytope type not supported.'));
end

% ------------------------------ END OF CODE ------------------------------
