function example_hybrid_reach_ARCH21_brake_BRKDC01()
% example_hybrid_reach_ARCH21_brake_BRKDC01 -  brake benchmark from the
%    2021 ARCH competition
%
% Syntax:  
%    example_hybrid_reach_ARCH21_brake_BRKDC01
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% 
% Author:       Niklas Kochdumper
% Written:      27-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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

% system matrices
A = [-(R+K^2/d_rot)/L, 0, K_P/L, K_I/L; ....
     K/(i*d_rot), 0, 0, 0; ...
     0, 0, 0, 0; ...
     0, 0, 0, 0];
 
B = [0;0;0;0];

% reset
reset.A = [1 0 0 0; 0 1 0 0; 0 -1 0 0;0 -T_sample 0 1];
reset.b = [0; 0; x0; T_sample*x0];

% linear system object
sys = linearSys(A,B);


% Parameter ---------------------------------------------------------------

params.tFinal = T_sample;
params.R0 = zonotope([0;0;0;0]);


% Settings ----------------------------------------------------------------

options.timeStep = T_sample/5;
options.taylorTerms = 5;
options.zonotopeOrder = 20;


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = {};

for i = 1:1000
    
   % compute reachable set
   Rtemp = reach(sys,params,options);
   R = [R; Rtemp.timeInterval.set];
   
   % compute reset
   params.R0 = reset.A*Rtemp.timePoint.set{end} + reset.b;
    
end

tComp = toc(clock);


% Verification ------------------------------------------------------------

clock = tic;
res = 1;

for i = 1:length(R)
   temp = supremum(interval(project(R{i},2)));
   if temp > x0
      res = 0;
      break;
   end
end
tVer = toc(clock);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tVer + tComp)]);
 

% Visualization -----------------------------------------------------------

figure; hold on

% plot reachable set
t = 0;

for i = 1:length(R)
    
    % get intervals
    intX = interval(project(R{i},2));
    intT = interval(t,t + options.timeStep);
    
    int = cartProd(intT,intX);
    
    % plot the interval
    plot(int,[1,2],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor',[.6 .6 .6]);
    
    % update time
    t = t + options.timeStep;
end

% plot specification
plot([0 0.1],[x0 x0],'--r');

% formatting
box on
xlim([0,0.1]);
ylim([0,0.06]);
xlabel('t');
ylabel('x');

%------------- END OF CODE --------------