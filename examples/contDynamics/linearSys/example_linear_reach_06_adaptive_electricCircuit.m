function res = example_linear_reach_06_adaptive_electricCircuit
% example_linear_reach_06_adaptive_electricCircuit - example for the
%    adaptive reachability algorithm from [1] using a simple 2D system
%
% Syntax:
%    res = example_linear_reach_06_adaptive_electricCircuit
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Wetzlinger et al. "Fully-automated verification of linear
%        systems using inner- and outer-approximations of reachable sets",
%        arXiv, 2022.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       22-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 0.02;
params.R0 = zonotope(interval([1;3],[3;5]));
params.U = zonotope(interval(-0.1,0.1));


% Reachability Settings ---------------------------------------------------

options.linAlg = 'adaptive';


% System Dynamics ---------------------------------------------------------

L = 2.5e-3;                         % inductance
C = 1.5e-3;                         % capacity
R = 2;                              % resistance

A = [-1/(R*C) 1/C; -1/L 0];
B = [0; 1/L];
dim_x = length(A);

sys = linearSys(A,B);


% Reachability Analysis ---------------------------------------------------

% with varying maximum error
errors = [0.04;0.02;0.01];
R = cell(length(errors),1);

for i = 1:length(errors)
    options.error = errors(i);
    clock = tic;
    R{i} = reach(sys,params,options);
    tComp = toc(clock);
    disp(['Computation time (',num2str(errors(i)),'): ',num2str(tComp),'s']);
end

fprintf('Compute inner-approximation... ');
% compute inner-approximation
pgon = cell(length(R),1);
for i = length(R):-1:1
    w = warning(); warning('off');
    V = vertices(R{i}.timePoint.set{end});
    pgon1 = polygon(V(1,:),V(2,:));
    V = sqrt(dim_x)*R{i}.timePoint.error(end)*[1 0 -1 0;0 -1 0 1];
    pgon2 = polygon(V(1,:),V(2,:));
    pgon{i,1} = minkDiff(pgon1,pgon2);
    warning(w);
end
fprintf('Done!\n');

% compute near-exact solution (very small error)
options.error = 0.0005;
Rexact = reach(sys,params,options);


% Visualization -----------------------------------------------------------

colors = {colorblind('b'),colorblind('r'),colorblind('g')};
han = []; text = cell(length(colors),1);

figure;
for p=1:2
    subplot(1,2,p); hold on; box on;
    % plot near-exact solution
    plot(Rexact.timePoint.set{end},[1,2],'k','LineWidth',1.2);
    for i=1:3
        % plot outer-approximation
        temp = plot(R{i}.timePoint.set{end},[1,2],...
            'EdgeColor',colors{i},'LineWidth',1.2);
        han = [han,temp];
        text{i} = ['$\varepsilon_{\max}$ = ',num2str(errors(i))];
        % plot inner-approximation
        plot(pgon{i},[1,2],'EdgeColor',colors{i},'LineWidth',1.2);
    end

    if p == 1
        axis([-0.4,0.15,-0.35,0.2]);
        legend(han,text{:},'interpreter','latex','Location','northwest');
        xlabel('$u_C$','interpreter','latex');
        ylabel('$i_L$','interpreter','latex');
    elseif p == 2
        % zoom and no ticks
        axis([0.05,0.15,-0.05,0.05]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
    end
end

% completed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
