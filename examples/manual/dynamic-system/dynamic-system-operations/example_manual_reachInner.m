function example_manual_reachInner()
% example_manual_reachInner - example from the manual demonstrating the 
% reachInner operation as defined in the manual
%
% Syntax:
%   example_manual_reachInner()
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

% system dynamics
f = @(x,u) [1-2*x(1) + 3/2*x(1)^2*x(2); ...
    x(1)-3/2*x(1)^2*x(2)];
sys = nonlinearSys(f);

% parameter
params.tFinal = 1;
params.R0 = interval([0.75;0],[1;0.25]);

% reachability settings
options.algInner = 'scale';
options.timeStep = 0.001;
options.taylorTerms = 10;
options.zonotopeOrder = 50;
options.intermediateOrder = 20;
options.errorOrder = 10;

% reachability analysis
[Rin,Rout] = reachInner(sys,params,options);

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(Rout.timePoint.set{end})
plot(Rin.timePoint.set{end})

enlargeAxis(1.2);
% title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')
xlim([0.6 0.8]); ylim([0.46 0.64])

% ------------------------------ END OF CODE ------------------------------
