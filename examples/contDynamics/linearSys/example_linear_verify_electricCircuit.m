function res = example_linear_verify_electricCircuit
% example_linear_verify_electricCircuit - application of the automated
%    verification algorithm [1] on an electric circuit [Sec. VII.-A, 1]
%
% Syntax:
%    res = example_linear_verify_electricCircuit
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Wetzlinger, N. Kochdumper, S. Bak, M. Althoff. "Fully Automated
%        Verification of Linear Systems Using Inner and Outer
%        Approximations of Reachable Sets", TAC 2023.
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

% system dynamics 
L = 2.5e-3;                         % inductance
C = 1.5e-3;                         % capacity
R = 2;                              % resistance

A = [-1/(R*C) 1/C; -1/L 0];
B = [0; 1/L];
dim_x = length(A);

sys = linearSys(A,B);

% reachability parameters
params.tFinal = 0.02;
params.R0 = zonotope(interval([1;3],[3;5]));
params.U = zonotope(interval(-0.1,0.1));
options.linAlg = 'adaptive';

% reachability analysis with varying maximum error
errors = {0.04;0.02;0.01};
R = cell(length(errors),1);

for i = 1:length(errors)
   options.error = errors{i};
   clock = tic;
   R{i} = reach(sys,params,options);
   tComp = toc(clock);
   disp(['Computation time (',num2str(errors{i}),'):',num2str(tComp),'s']);
end

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

% compute near-exact solution (very small error)
options.error = 0.0005;
clock = tic;
Rexact = reach(sys,params,options);
tComp = toc(clock);

% visualization
colors = {colorblind('b'),colorblind('r'),[0 0.8 0]};
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
        han = [temp,han];
        text{i} = ['$\varepsilon_{\max}$ = ',num2str(errors{i})];
        % plot inner-approximation
        plot(pgon{i},[1,2],'EdgeColor',colors{i},'LineWidth',1.2);
    end

    if p == 1
        axis([-0.4,0.15,-0.35,0.2]);
        legend(han,text{:},'interpreter','latex',...
            'Location','northwest','FontSize',14);
        xlabel('$u_C$','interpreter','latex','FontSize',16);
        ylabel('$i_L$','interpreter','latex','FontSize',16);
    end
    if p == 2
        % zoom and no ticks
        axis([0.05,0.15,-0.05,0.05]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
    end
end

% completed
res = true;

% ------------------------------ END OF CODE ------------------------------
