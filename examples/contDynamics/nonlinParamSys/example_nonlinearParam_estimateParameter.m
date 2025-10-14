function completed = example_nonlinearParam_estimateParameter()
% example_nonlinearParam_estimateParameter - example of paramter estimation 
%    for nonlinearParam systems.
%
% Syntax:
%    completed = example_nonlinearParam_estimateParameter()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Authors:       Laura Luetzow
% Written:       25-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
p_true = 2.8;
f = @(x,u,p) [x(1)*cos(x(2)); ...
    x(1)^2/p(1) * tan(u(1))];
sys = nonlinParamSys(f);

% simulate ground-truth system
simOpts.x0 = [1; 0];
simOpts.tFinal = 2;
simOpts.u = 0.3;
simOpts.p = p_true;
[t,x] = simulate(sys,simOpts);
u = repmat(simOpts.u,[1,length(t)]);

p_est = estimateParameter(sys,x,t,u);

% simulate system with estimated parameter
simOpts.p = p_est;
[t_est,x_est] = simulate(sys,simOpts);

% Visualization 

% plot different projections
for i = 1:2
    
    figure; hold on; box on;
    useCORAcolors("CORA:contDynamics", 2)

    % plot reachable sets
    plot(t,x(i,:),'DisplayName','true System');
    plot(t_est,x_est(i,:),'r','DisplayName','estimated');
    
    % label plot
    xlabel('t');
    ylabel(['x_{',num2str(i),'}']);
    legend();
end

completed = true;

% ------------------------------ END OF CODE ------------------------------
