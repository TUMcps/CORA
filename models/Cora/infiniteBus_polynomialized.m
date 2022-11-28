function f = infiniteBus_polynomialized(x,u)
% infiniteBus - infinite-bus power system from [1] polynomialized according
% to [2].
%
% Syntax:  
%    f = infiniteBus(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%    [1] M. Althoff. "Formal Verification of Power Systems", submitted to
%        ARCH 2022.
%    [2] M. Anghel, F. Milano and A. Papachristodoulou, "Algorithmic 
%        Construction of Lyapunov Functions for Power System Stability 
%        Analysis," in IEEE Transactions on Circuits and Systems I: Regular 
%        Papers, vol. 60, no. 9, pp. 2533-2546, Sept. 2013, 
%        doi: 10.1109/TCSI.2013.2246233.

% Author:       Matthias Althoff
% Written:      02-Jun-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% dynamics is eqiuivalent to "infiniteBus":
% \dot{x}_1 = x_2;
% \dot{x}_2 = 15*pi*(1 - sin(x_1) - 0.04*x_2);
%
% change of variables:
% z_1 = sin(x_1)
% z_2 = 1 - cos(x_1)
% z_3 = x_2
%
% new dynamics:
% \dot{z}_1 = z_3 - z_2*z_3;
% \dot{z}_2 = z_1*z_3;
% \dot{z}_3 = 15*pi*(1 - z_1 - 0.04*z_3);

    f = [x(3) - x(2)*x(3);...
         x(1)*x(3); ...
         15*pi*(1 - x(1) - 0.04*x(3))];
     
     
end

%------------- END OF CODE --------------