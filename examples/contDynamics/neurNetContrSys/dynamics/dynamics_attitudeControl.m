function f = dynamics_attitudeControl(x,u)
% dynamics_attitudeControl - attitude control of a rigid body [1]
%
% Syntax:
%    f = dynamics_attitudeControl(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
% 
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] Prajna et al.: Nonlinear control synthesis by convex optimization. 
%       IEEE Transactions on Automatic Control 49(2), 310–314 (2004)

% Authors:       Tobias Ladner
% Written:       20-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parameters

    % extract variables for better readability
    w1 = x(1);
    w2 = x(2);
    w3 = x(3);
    p1 = x(4);
    p2 = x(5);
    p3 = x(6);

    u0 = u(1);
    u1 = u(2);
    u2 = u(3);

    % dynamic equations
    
    f = [
         0.25*(u0 + w2*w3);
         0.5*(u1-3*w1*w3);
         u2 + 2*w1*w2;
         0.5*(w2*(p1^2+p2^2+p3^2-p3)+w3*(p1^2+p2^2+p2+p3^2)+w1*(p1^2+p2^2+p3^2+1))
         0.5*(w1*(p1^2+p2^2+p3^2+p3)+w3*(p1^2-p1+p2^2+p3^2)+w2*(p1^2+p2^2+p3^2+1))
         0.5*(w1*(p1^2+p2^2-p2+p3^2)+w2*(p1^2+p1+p2^2+p3^2)+w3*(p1^2+p2^2+p3^2+1))
    ];
     
end

% ------------------------------ END OF CODE ------------------------------
