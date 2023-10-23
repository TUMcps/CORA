function text = benchmark_hybrid_reach_ARCH23_brake_BRKNP01()
% benchmark_hybrid_reach_ARCH23_brake_BRKNP01 -  brake benchmark from the
%    2023 ARCH competition (parametric + nondet. switching)
%
% Syntax:
%    benchmark_hybrid_reach_ARCH23_brake_BRKNP01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%

% Authors:       Niklas Kochdumper
% Written:       27-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Dynamic System ----------------------------------------------------------

% parameter values
L = 0.001;
K_P = 10000;
K_I = 1000;
R = 0.5;
K = 0.02;
d_rot = 0.1;
i = 113.1167;
T_sample = 0.0001;
x0 = 0.05;
t_min = -1e-8;
t_max = 1e-7;

% system matrices
Ac = [-(R+K^2/d_rot)/L, 0, K_P/L, K_I/L, 0; ....
      K/(i*d_rot), 0, 0, 0, 0; ...
      0, 0, 0, 0, 0; ...
      0, 0, 0, 0, 0; ...
      0, 0, 0, 0, 0];
Ag{1} = [3 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
A = matZonotope(Ac,Ag);
 
B = [0;0;0;0;1];

% reset
reset.A = [1 0 0 0 0; 0 1 0 0 0; 0 -1 0 0 0;0 -T_sample 0 1 0;0 0 0 0 1];
reset.c = [0; 0; x0; T_sample*x0; -T_sample];

% invariant set
inv = polytope([0 0 0 0 1],T_sample + t_max);

% guard set
guard = polytope([0 0 0 0 1; 0 0 0 0 -1], ...
                    [T_sample+t_max;-(T_sample+t_min)]);

% linear paramtric system
sys = linParamSys(A,B,'constParam');

% hybrid automaton
trans(1) = transition(guard,reset,1);
loc(1) = location(inv,trans,sys);

HA = hybridAutomaton(loc);


% Parameter ---------------------------------------------------------------

params.tFinal = 0.1;
params.R0 = zonotope([0;0;0;0;0],diag(1e-8*[1 1 1 1 1]));
params.U = zonotope(1);
params.startLoc = 1;


% Settings ----------------------------------------------------------------

options.timeStep = T_sample/5;
options.taylorTerms = 10;
% options.intermediateOrder = 3;
options.zonotopeOrder = 20;
options.intermediateTerms = 3;

options.guardIntersect = 'conZonotope';
options.enclose = {'box','pca'};
options.guardOrder = 5;


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = reach(HA,params,options);
tComp = toc(clock);


% Verification ------------------------------------------------------------

tMax = 0;
t = 0;

for i = 1:size(R,1)
    for j = 1:length(R(i).timeInterval.set)
       
       % get intervals
       temp = R(i).timeInterval.set{j};
       intX = interval(project(temp,2));
       intT = interval(project(temp,5)) + t;
        
       % check specification
       if supremum(intX) < x0
          tMax = supremum(intT);
          break;
       end
    end
    
    % update time
    t = t + T_sample;
end

res = tMax < 0.1;

disp(['max. time horizon: ',num2str(tMax)]);
disp(['computation time: ', num2str(tComp)]);
 

% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")

% plot reachable set
plotOverTime(R,2);

% plot specification
spec = specification(interval([0; x0], [0.1; x0]));
plot(spec,[1 2],'--');

% formatting
xlim([0,0.1]);
ylim([0,0.06]);
xlabel('t');
ylabel('x');

text = ['Brake,BRKNP01,',num2str(res),',',num2str(tComp)];

% ------------------------------ END OF CODE ------------------------------
