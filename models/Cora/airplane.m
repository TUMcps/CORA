function f = airplane(x,u)
% airplane - airplane from Section 3.7 in [1]
%
% Syntax:  
%    f = airplane(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] T. Johnson and et al. "ARCH-COMP21 Category Report: Artificial 
%       Intelligence and Neural Network Control Systems (AINNCS) for 
%       Continuous and Hybrid Systems Plants", 2021

% Author:       Niklas Kochdumper
% Written:      08-December-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parameter
    m = 1; 
    g = 1;
    I_x = 1;
    I_y = 1;
    I_z = 1;
    I_xz = 0;

    % dynamic equations (states)
    T1 = [cos(x(9)) -sin(x(9)) 0; sin(x(9)) cos(x(9)) 0; 0 0 1];
    T2 = [cos(x(8)) 0 sin(x(8)); 0 1 0; -sin(x(8)) 0 cos(x(8))];
    T3 = [1 0 0; 0 cos(x(7)) -sin(x(7)); 0 sin(x(7)) cos(x(7))];
    
    f1 = T1*T2*T3*[x(4);x(5);x(6)];
    
    % dynamic equations (velocities)
    f2 = [-g*sin(x(8)) + u(1)/m - x(12)*x(6) + x(10)*x(5); ...
          g*cos(x(8))*sin(x(7)) + u(2)/m - x(10)*x(4) + x(11)*x(6); ...
          g*cos(x(8))*cos(x(7)) + u(3)/m - x(11)*x(5) + x(12)*x(4)];
      
    % dynamic equations (angles)
    f3 = [1 tan(x(8))*sin(x(7)) tan(x(8))*cos(x(7)); 
          0 cos(x(7)) -sin(x(7));
          0 sin(x(7))/cos(x(8)) cos(x(7))/cos(x(8))]*[x(11);x(12);x(10)];
      
    % dynamic equations (rotation rates)
    a = 1/(I_xz^2 - I_x*I_z);
    M = [I_xz/a -I_x/a 0;-I_z/a I_xz/a 0; 0 0 I_y];
    
    f4 = M * [u(4) - (I_z-I_y)*x(12)*x(10) - I_xz*x(11)*x(12); ...
              u(6) - (I_y-I_x)*x(12)*x(11) - I_xz*x(10)*x(12); ...
              u(5) - I_xz*(x(10)^2-x(11)^2) - (I_x-I_z)*x(11)*x(10)];
    
    % overall dynamic equations
    f = [f1; f2; f3; f4];
     
end

%------------- END OF CODE --------------