function res = example_nonlinear_reach_10_pendulum()
% example_nonlinear_reach_10_pendulum - example for nonlinear reachability
%    analysis: a simple pendulum swinging in a 2D plane
%
% Syntax:
%    example_nonlinear_reach_10_pendulum
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       03-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tStart = 0;   % Start time
params.tFinal = 2;   % Final time

P.L = 0.2;           % length of the rod (m)
P.mu = 0.5;          % friction


% Initial set for reachability analysis
initangle = 120; % degrees
initangle = initangle * (2 * pi / 360); % convert to radians
initangledelta = 2; % degrees
initangledelta = initangledelta * (2 * pi / 360); % convert to radians
initvelocity = 0; % m/s
params.R0 = zonotope([initangle;initvelocity],[initangledelta,0;0,0]);

% Input set
params.U = zonotope(0); % no external inputs


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.005;        % Time step size
options.taylorTerms = 5;        % Taylor terms
options.zonotopeOrder = 30;     % Zonotope order

options.alg                  = 'lin';
options.tensorOrder          = 3;
options.errorOrder           = 5;
options.intermediateOrder    = 10;
options.maxError             = Inf(dim(params.R0),1);
options.reductionInterval    = Inf;


% System Dynamics ---------------------------------------------------------

dynFun = @(x,u) pendulum(x,u,P);
sys = nonlinearSys('pendulum',dynFun,2,1);


% Reachability Analysis ---------------------------------------------------

% Compute reachable set 
R = reach(sys, params, options);


% Simulation --------------------------------------------------------------

simOpt.points = 20;
simOpt.fracVert = 0.1;
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

% use pause() function to mimic video

% set up figure (label, title, axis)
figure;
subplot(1,2,1); hold on; box on;
xlabel("Angle"); ylabel("Angular Velocity"); title("Phase Space");
axis([-pi,pi,-15,15]);
subplot(1,2,2); hold on; box on;
xlabel("x"); ylabel("y"); title("XY Space");
axis([-P.L*1.5,P.L*1.5,-P.L*1.5,P.L*1.5]);
% plot hinge of pendulum
scatter(0,0,16,'r','filled');

% convert phase space to x-y plane
% pendulum locked on origin, length L
% ... x = L * sin phi, y = -L * cos phi
% note: usage of intervals does not contain restriction to circle equation!
Rphasespace = query(R,'reachSet');
totalSets = length(Rphasespace);
xyInt = cell(totalSets,1);
% reachable sets
for i=1:length(Rphasespace)
    temp = interval(Rphasespace{i});
    angleminmax = temp(1);
    xyInt{i,1} = [P.L*sin(angleminmax);-P.L*cos(angleminmax)];
end
% simulation
simX = cell(simOpt.points,1);
simY = cell(simOpt.points,1);
idx = cell(simOpt.points,1);
for r=1:simOpt.points
    phi = simRes(r).x{1}(:,1);
    simX{r,1} = P.L*sin(phi);
    simY{r,1} = -P.L*cos(phi);
    % prepare indices for simulation
    idx{r} = [1;zeros(totalSets,1)];
    for i=1:totalSets
        idx{r}(i+1) = find(simRes(r).t{1} <= options.timeStep*i,1,'last');
    end
end

for i=1:totalSets
    % grayscale color brightening over time
    grayscale = CORAcolor("CORA:reachSet") * i/totalSets;
    
    % phase space plot
    subplot(1,2,1);
    % plot reachable set
    plot(Rphasespace{i},[1,2],'FaceColor',grayscale);
    % plot simulation
    for r=1:simOpt.points
        % always plot until last index before current time
        % add idx below!
        plot(simRes(r).x{1}(idx{r}(i):idx{r}(i+1),1),...
            simRes(r).x{1}(idx{r}(i):idx{r}(i+1),2),'Color',CORAcolor("CORA:simulations"));
    end
    
    % xy space plot (additional over-approximation in conversion)
    subplot(1,2,2);
    % reachable set
    plot(xyInt{i},[1,2],'FaceColor',grayscale);
    % simulation
    for r=1:simOpt.points
        % take same indices as above
        plot(simX{r}(idx{r}(i):idx{r}(i+1)),...
            simY{r}(idx{r}(i):idx{r}(i+1)),'Color',CORAcolor("CORA:simulations"));
    end
    
    % pause for video effect
    pause(0.2);
end


% example completed
res = true;   

end

% ------------------------------ END OF CODE ------------------------------
