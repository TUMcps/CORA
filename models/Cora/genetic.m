function dx = genetic(x,u)
% genetic - genetic model from [2,V.]
% 
%           name  | init set [1]  | init set [2]
%   x(1):   DA      1               1
%   x(2):   DR      1               1.3
%   x(3):   D'A     0               0.1
%   x(4):   D'R     0               0.1
%   x(5):   MA      0               0.1
%   x(6):   A       0               1.3
%   x(7):   MR      0               2.5
%   x(8):   R       0               0.6
%   x(9):   C       0               1.3
%                   (bio.paper)     -> all -+ W0 = 0.01|0.02|0.04
%
% Syntax:  
%    dx = genetic(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
%
% Reference:
%   [1] J. Vilar, H. Kueh, N. Barkai, S. Leibler. "Mechanisms of
%           noise-resistance in genetic oscillators".
%   [2] X. Chen, S. Sankaranarayanan. "Decomposed Reachability Analysis
%           for Nonlinear Systems".

% Author:        Mark Wetzlinger
% Written:       10-July-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
 
% constants     unit                [1]     [2]
thetaA = 50;    % h^-1              50      50
gammaA = 1;     % mol^-1 h^-1       1       0.1
thetaR = 100;   % h^-1              100     100
gammaR = 1;     % mol^-1 h^-1       1       1
alpha_A = 500;  % h^-1              500     5
alphaA = 50;    % h^-1              50      0.5
deltaMA = 10;   % h^-1              10      10
betaA = 50;     % h^-1              50      50
gammaC = 2;     % mol^-1 h^-1       2       2
deltaA = 1;     % h^-1              1       1
alpha_R = 50;   % h^-1              50      50
alphaR = 0.01;  % h^-1              0.01    0.01
deltaMR = 0.5;  % h^-1              0.5     0.5
betaR = 5;      % h^-1              5       0.5
deltaR = 0.2;   % h^-1              0.2     0.2

dx(1,1) = thetaA*x(3) - gammaA*x(1)*x(6);
dx(2,1) = thetaR*x(4) - gammaR*x(2)*x(6);
dx(3,1) = gammaA*x(1)*x(6) - thetaA*x(3);
dx(4,1) = gammaR*x(2)*x(6) - thetaR*x(4);
dx(5,1) = alpha_A*x(3) + alphaA*x(1) - deltaMA*x(5);
dx(6,1) = betaA*x(5) + thetaA*x(3) + thetaR*x(4) - ...
    x(6)*(gammaA*x(1) + gammaR*x(2) + gammaC*x(8) + deltaA);
dx(7,1) = alpha_R*x(4) + alphaR*x(2) - deltaMR*x(7);
dx(8,1) = betaR*x(7) - gammaC*x(6)*x(8) + deltaA*x(9) - deltaR*x(8);
dx(9,1) = gammaC*x(6)*x(8) - deltaA*x(9);
    
end

%------------- END OF CODE --------------