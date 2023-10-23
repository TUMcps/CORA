function example_manual_simulate()
% example_manual_simulate - example from the manual demonstrating the 
% simulate operation as defined in the manual
%
% Syntax:
%   example_manual_simulate()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% nonlinear system
f = @(x,u) [x(2) + u; ...
    (1-x(1)^2)*x(2)-x(1)];
sys = nonlinearSys(f);

% parameter
params.x0 = [1.4;2.3];
params.tFinal = 6;
params.u = [0.1 0 -0.1 0.2];

% simulation
[t,x] = simulate(sys,params);

% plot --------------------------------------------------------------------

figure; hold on;
plot(x(:,1),x(:,2),'k')
plot(x(1,1),x(1,2),'.r')

% title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
