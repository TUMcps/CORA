function text = benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISU01
% benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISU01 - iss benchmark from
%     the 2023 ARCH competition
%
% Syntax:
%    text = benchmark_linear_verifyFast_ARCH23_iss_ISSF01_ISU01()
%
% Inputs:
%    -
%
% Outputs:
%    text - char array

% Authors:       Mark Wetzlinger
% Written:       23-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

% initial set
R0 = interval(-0.0001*ones(270,1),0.0001*ones(270,1));
params.R0 = zonotope(R0);

% uncertain inputs
U = interval([0;0.8;0.9],[0.1;1;1]);
params.U = zonotope(U);

% final time
params.tFinal = 20;

options = struct();
options.verifyAlg = 'reachavoid:supportFunc';

% Specification -----------------------------------------------------------

% forall t: -5e-4 <= y3 <= 5e-4
d = 5e-4;
hs1 = halfspace([0, 0, 1], -d);
hs2 = halfspace([0, 0, -1], -d);
spec = specification({hs1,hs2},'unsafeSet');


% System Dynamics ---------------------------------------------------------

% load system matrices
load iss.mat A B C

% construct the linear system object
sys = linearSys('iss',A,B,[],C);


% Verification ------------------------------------------------------------

% min steps needed: 150
[res,fals,savedata] = verify(sys,params,options,spec);


% Return value ------------------------------------------------------------

text = ['Spacestation,ISSF01_ISU01,',num2str(res),',',num2str(savedata.tComp)];

% Visualization -----------------------------------------------------------

% figure;
% 
% % convert to discrete-time system for faster evaluation
% sysDT = linearSysDT(sys,fals.tu(2));
% [t,~,~,y] = simulate(sysDT,fals);
% 
% % project initial set and trajectory on the halfspace of the specification
% Y0_proj = interval( (hs1.c' * full(sys.C)) * params.R0 );
% Y0_proj_inf = infimum(Y0_proj);
% Y0_proj_sup = supremum(Y0_proj);
% y_proj = (hs1.c' * y')';
% % plot over time
% subplot(2,1,1); hold on; box on;
% xlabel('$t$','interpreter','latex');
% ylabel('$y_3$','interpreter','latex');
% xlim([0,15]);
% ylim([-6e-4,6e-4]);
% % falsifying trajectory
% h_y = plot(t,y_proj,'Color',colorblind('b'));
% % spec
% h_spec = plot([0,params.tFinal],[d,d],'Color',colorblind('r'),'LineStyle','--');
% h_spec = plot([0,params.tFinal],[-d,-d],'Color',colorblind('r'),'LineStyle','--');
% 
% % figure showing input trajectory
% subplot(2,1,2); hold on; box on;
% xlabel('$t$','interpreter','latex');
% ylabel('$u_i(t)$','interpreter','latex');
% ylim([-0.1,1.1]);
% xlim([0,15]);
% 
% deltat = fals.tu(2) - fals.tu(1);
% tVec{1} = [fals.tu(1);repelem(fals.tu(2:end),2);fals.tu(end)+deltat];
% tVec{2} = [fals.tu(1);repelem(fals.tu(2:end),2);fals.tu(end)+deltat];
% tVec{3} = [fals.tu(1);repelem(fals.tu(2:end),2);fals.tu(end)+deltat];
% u{1} = repelem(fals.u(1,:),2)';
% u{2} = repelem(fals.u(2,:),2)';
% u{3} = repelem(fals.u(3,:),2)';
% % remove duplicates
% for uu=1:3
%     idx{uu} = true(length(u{uu}),1);
%     i = 2;
%     while i < length(u{uu})
%         % check if next index has the same value
%         if withinTol(u{uu}(i-1),u{uu}(i),eps) && withinTol(u{uu}(i),u{uu}(i+1),eps)
%             idx{uu}(i) = false;
%         end
%         i = i+1;
%     end
%     tVec{uu} = tVec{uu}(idx{uu});
%     u{uu} = u{uu}(idx{uu});
% end
% 
% h_u1 = plot(tVec{1},u{1},'Color',[0 0.5 0]);
% h_u2 = plot(tVec{2},u{2},'Color',colorblind('r'));
% h_u3 = plot(tVec{3},u{3},'k');
% 
% legend([h_u1,h_u2,h_u3],'$u_1$','$u_2$','$u_3$','interpreter','latex',...
%     'Orientation','horizontal','Location','east');


% ------------------------------ END OF CODE ------------------------------
