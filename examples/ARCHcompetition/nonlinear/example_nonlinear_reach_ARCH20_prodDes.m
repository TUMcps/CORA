function example_nonlinear_reach_ARCH20_prodDes()
% example_nonlinear_reach_ARCH20_prodDes - production-destruction benchmark
%                                          from the 2020 ARCH competition
%
% Syntax:  
%    example_nonlinear_reach_ARCH20_prodDes()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      17-May-2020
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


% Parameter ---------------------------------------------------------------

params.tFinal = 100;

% initial set
R0_I = zonotope(interval([9.5;0.01;0.01],[10;0.01;0.01]));
R0_P = zonotope(interval([9.98;0.01;0.01;0.296],[9.98;0.01;0.01;0.304]));
R0_IP = zonotope(interval([9.7;0.01;0.01;0.298],[10;0.01;0.01;0.302]));


% Reachability Options ----------------------------------------------------

% general options
options.zonotopeOrder = 20;
options.taylorTerms = 10;

options.errorOrder = 5;
options.intermediateOrder = 10;

options.lagrangeRem.simplify = 'optimize';

% select algorithm
options.alg = 'lin';
options.tensorOrder = 3;


% System Dynamics ---------------------------------------------------------

sys_I = nonlinearSys(@prodDes);
sys_P = nonlinearSys(@prodDesParam);
sys_IP = sys_P;


% Reachability Analysis ---------------------------------------------------

% reachability analyis (case I)
params.R0 = R0_I;
options.timeStep = params.tFinal/1700;

tic;
R_I = reach(sys_I,params,options);
tComp_I = toc;

% reachability analyis (case P)
params.R0 = R0_P;
options.timeStep = params.tFinal/3000;

tic;
R_P = reach(sys_P,params,options);
tComp_P = toc;

% reachability analyis (case IP)
params.R0 = R0_IP;
options.timeStep = params.tFinal/3000;

tic;
R_IP = reach(sys_IP,params,options);
tComp_IP = toc;



% Verification ------------------------------------------------------------

% verification (case I)
res1 = 1; res2 = 1;

int = interval(R_I.timePoint.set{end});
infi = infimum(int);

if infi(1) < 0 || infi(2) < 0 || infi(3) < 0
    res1 = 0; 
end

Rtemp = [1 1 1]*R_I.timePoint.set{end};
intSum = interval(Rtemp);

if ~in(intSum,10)
    res2 = 0; 
end

disp(' ');
disp('Case I ------------------')
disp(['Specifications satisfied: ',num2str(res1 & res2)]);
disp(['Computation time: ',num2str(tComp_I),'s']);
disp(['Volume: ',num2str(volume(int))]);
disp(' ');

% verification (case P)
res1 = 1; res2 = 1;

int = interval(R_P.timePoint.set{end});
infi = infimum(int);

if infi(1) < 0 || infi(2) < 0 || infi(3) < 0
    res1 = 0; 
end

Rtemp = [1 1 1 0]*R_P.timePoint.set{end};
intSum = interval(Rtemp);

if ~in(intSum,10)
    res2 = 0; 
end

disp('Case P ------------------')
disp(['Specifications satisfied: ',num2str(res1 & res2)]);
disp(['Computation time: ',num2str(tComp_P),'s']);
disp(['Volume: ',num2str(volume(int))]);
disp(' ');

% verification (case IP)
res1 = 1; res2 = 1;

int = interval(R_IP.timePoint.set{end});
infi = infimum(int);

if infi(1) < 0 || infi(2) < 0 || infi(3) < 0
    res1 = 0; 
end

Rtemp = [1 1 1 0]*R_IP.timePoint.set{end};
intSum = interval(Rtemp);

if ~in(intSum,10)
    res2 = 0; 
end

disp('Case I&P ----------------')
disp(['Specifications satisfied: ',num2str(res1 & res2)]);
disp(['Computation time: ',num2str(tComp_IP),'s']);
disp(['Volume: ',num2str(volume(int))]);
disp(' ');



% Visualization -----------------------------------------------------------

figure; hold on; box on;

% visualization (case I)
plotOverTime(R_I,3,'b','EdgeColor','b');

% visualization (case IP)
plotOverTime(R_IP,3,'FaceColor',[0 .6 0],'EdgeColor',[0 .6 0]);

% visualization (case P)
plotOverTime(R_P,3,'r','EdgeColor','r');

xlabel('t')
ylabel('z')
xlim([0,100]);
ylim([0,11]);