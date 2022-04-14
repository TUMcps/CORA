function dx = tank6paramEq(x,u,p)
% tank6paramEq - system dynamics for the tank benchmark 
%                (see Sec. VII in [1])
%
% Syntax:  
%    dx = tank6paramEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%    p - parameter vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Author:        Matthias Althoff
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % parameter
    k = 0.01;
    g = 9.81; 

    % differential equation
    dx(1,1)=u(1)+0.1+k*(4-x(6))-p(1)*sqrt(2*g)*sqrt(x(1));  % tank 1
    dx(2,1)=p(1)*sqrt(2*g)*(sqrt(x(1))-sqrt(x(2)));         % tank 2
    dx(3,1)=p(1)*sqrt(2*g)*(sqrt(x(2))-sqrt(x(3)));         % tank 3
    dx(4,1)=p(1)*sqrt(2*g)*(sqrt(x(3))-sqrt(x(4)));         % tank 4
    dx(5,1)=p(1)*sqrt(2*g)*(sqrt(x(4))-sqrt(x(5)));         % tank 5
    dx(6,1)=p(1)*sqrt(2*g)*(sqrt(x(5))-sqrt(x(6)));         % tank 6
    
%------------- END OF CODE --------------