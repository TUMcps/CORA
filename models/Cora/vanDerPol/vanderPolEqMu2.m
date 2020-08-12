function dx = vanderPolEqMu2(x,u)
% vanDerPolEqMu2 - system dynamics for the Van-der-Pol oscillator with
%                  increased stiffness (see Sec. 3.1 in [1])
%
% Syntax:  
%    dx = vanDerPolEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] F. Immler, “ARCH-COMP19 Category Report: Continuous and Hybrid 
%        Systems with Nonlinear Dynamics", 2019

% Author:        Matthias Althoff
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    mu=2;

    dx(1,1)=x(2);
    dx(2,1)=mu*(1-x(1)^2)*x(2)-x(1)+u(1);
    
end

%------------- END OF CODE --------------