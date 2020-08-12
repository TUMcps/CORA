function f = cstrDiscr(x,u,T)
% cstrDiscr - discrete-time version of the stirred-tank reactor system in
%             Sec. 6 in [1]
%
% Syntax:  
%     f = cstrDiscr(x,u,T)
%
% Inputs:
%    x - state vector
%    u - input vector
%    T - sampling time
%
% Outputs:
%    f - new state vector
%
% Reference:
%   [1] J. M. Bravo. "Robust MPC of constrained discrete-time nonlinear 
%       systems based on approximated reachable sets", Automatica 2006

% Author:       Niklas Kochdumper
% Written:      19-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % Parameter
    rho = 1000;
    Cp = 0.239;
    deltaH = -5e4;
    E_R = 8750;
    k0 = 7.2e10;
    UA = 5e4;
    q = 100;
    Tf = 350;
    V = 100;
    C_Af = 1;
    C_A0 = 0.5;
    T_0 = 350;
    T_c0 = 300;

    % Control law
    U = [-3 -6.9] * x;

    % State offset
    x = x + [C_A0;T_0];

    % Control law offset
    U = U + T_c0;

    % Dynamics
    f(1,1) = ((1-(q*T)/(2*V) - k0*T*exp(-E_R/x(2)))*x(1) + q/V * C_Af * T)/ ...
             (1 + (q*T)/(2*V)) + u(1)*T;

    f(2,1) = (x(2)*(1-0.5*T- (T*UA)/(2*V*rho*Cp)) + T*(Tf*q/V + (UA*U)/(V*rho*Cp)) ...
              - x(1)*(deltaH*k0*T)/(rho*Cp) * exp(-E_R/x(2))) ...
              / (1+0.5*T*q/V+(T*UA)/(2*V*rho*Cp)) + u(2)*T;

    % Offset substraction
    f = f - [C_A0;T_0];

end

%------------- END OF CODE --------------