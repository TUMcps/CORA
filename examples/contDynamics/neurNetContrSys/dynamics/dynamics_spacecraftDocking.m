function f = dynamics_spacecraftDocking(x,u)
% dynamics_spacecraftDocking - neural-network controlled space craft docking [1]
%
% Syntax:
%    f = dynamics_spacecraftDocking(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
% 
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] Ravaioli, U.: Spacecraft Docking Benchmark (2022)

% Authors:       Tobias Ladner
% Written:       20-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parameter
m = 12;
n = 0.001027;

% f(x, u) = A*x + B*u
A = [0 0 1 0; 
     0 0 0 1; 
     3*n^2 0 0 2*n; 
     0 0 -2*n 0];
B = [0 0; 
     0 0; 
     1/m 0; 
     0 1/m;];

f = A*x(1:4,:) + B*u(1:2,:);

end

% ------------------------------ END OF CODE ------------------------------
