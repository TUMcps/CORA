function res = test_hybridAutomaton_simulate
% test_hybridAutomaton_simulate - test function for simulate
%
% Syntax:
%    res = test_hybridAutomaton_simulate
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate simple automaton:
% - 1st median is guard set
% - right part moves clockwise
% - left part moves counter-clockwise
inv1 = polytope([-1 1]/sqrt(2),0);
inv2 = polytope([1 -1]/sqrt(2),0);

guard = conHyperplane([-1 1]/sqrt(2),0);
reset1 = struct('A',eye(2),'c',[-1;1]);
trans1 = transition(guard,reset1,2);

reset2 = struct('A',eye(2),'c',[1;-1]);
trans2 = transition(guard,reset2,1);

% clockwise motion
dyn1 = linearSys([0 1; -1 0],1);
% counter-clockwise motion
dyn2 = linearSys([0 -1; 1 0],1);

loc1 = location('clockwise',inv1,trans1,dyn1);
loc2 = location('counter-clockwise',inv2,trans2,dyn2);
HA = hybridAutomaton([loc1;loc2]);

% model parameters
params.x0 = [2;-2];
params.startLoc = 1;
params.finalLoc = 0;
params.tFinal = 3;

% simulate trajectory
[t,x,loc] = simulate(HA,params); 

% must be cell-arrays
res = iscell(t) && iscell(x) && isnumeric(loc);
% must be of same length
res(end+1,1) = length(t) == length(x) && length(t) == length(loc);
% check whether points are contained in respective invariant
res(end+1,1) = all(contains_(inv1,vertcat(x{loc==1})','exact',1e-10));
res(end+1,1) = all(contains_(inv2,vertcat(x{loc==2})','exact',1e-10));
% time before and after jumps must be the same, but state not
for i=1:length(t)-1
    res(end+1,1) = withinTol(t{i}(end),t{i+1}(1));
    res(end+1,1) = ~compareMatrices(x{i}(end,:)',x{i+1}(1,:)');
end


% automaton with different number of states per location
inv1 = polytope([-1 1]/sqrt(2),0);
inv2 = polytope([1 -1 0]/sqrt(2),0);

guard1 = conHyperplane([-1 1]/sqrt(2),0);
guard2 = conHyperplane([-1 1 0]/sqrt(2),0);
reset1 = struct('A',[1 0; 0 1; 0 0],'c',[-1;1;1]);
trans1 = transition(guard1,reset1,2);

reset2 = struct('A',[1 0 0; 0 1 0],'c',[1;-1]);
trans2 = transition(guard2,reset2,1);

% clockwise motion
dyn1 = linearSys([0 1; -1 0],1);
% counter-clockwise motion
dyn2 = linearSys([0 -1 0; 1 0 0; 0 0 0],1);

loc1 = location('clockwise',inv1,trans1,dyn1);
loc2 = location('counter-clockwise',inv2,trans2,dyn2);
HA = hybridAutomaton([loc1;loc2]);

% model parameters
params.x0 = [2;-2];
params.startLoc = 1;
params.finalLoc = 0;
params.tFinal = 3;

% simulate trajectory
[t,x,loc] = simulate(HA,params); 

% must be cell-arrays
res(end+1,1) = iscell(t) && iscell(x) && isnumeric(loc);
% must be of same length
res(end+1,1) = length(t) == length(x) && length(t) == length(loc);
% check whether points are contained in respective invariant
res(end+1,1) = all(contains_(inv1,vertcat(x{loc==1})','exact',1e-10));
res(end+1,1) = all(contains_(inv2,vertcat(x{loc==2})','exact',1e-10));
% time before and after jumps must be the same, but state not
for i=1:length(t)-1
    res(end+1,1) = withinTol(t{i}(end),t{i+1}(1));
    res(end+1,1) = ~compareMatrices(x{i}(end,:)',x{i+1}(1,:)');
end


% automaton with nonlinear invariant and guard sets (previously there was
% an error for invariants with multiple inequality constraints, which is
% the case we are testing here)
A = [0.01 0;0 0.01];
B = [0;0];
c = [-1; 0];
linSys = linearSys(A,B,c);

x = sym('x',[2,1]);
inv = levelSet([x(1) - 1; x(1).^2 + x(2).^2 - 4],x,'<=');

guard = levelSet(x(1) - 1,x,'==');
reset.A = zeros(2); reset.c = zeros(2,1);
trans = transition(guard,reset,1);

guard = levelSet(x(1).^2 + x(2).^2 - 4,x,'==');
reset.A = zeros(2); reset.c = zeros(2,1);
trans = [trans;transition(guard,reset,1)];

HA = hybridAutomaton(location(inv,trans,linSys));

% simulation
params.x0 = [0;0];
params.startLoc = 1;
params.finalLoc = 0;
params.tFinal = 4;

[~,x,~] = simulate(HA,params); 

% check if trajectory stayed within the invariant
for i = 1:length(x)
    for j = 1:size(x{i},1)
        if ~contains(inv,x{i}(j,:)')
            res(end+1,1) = false;
        end
    end
end


% test the case with time-varying inputs by comparing the simulation of a
% continuous system with the output of an equivalent automaton
linSys = linearSys([-0.7 -2; 2 -0.7],1);

R0 = zonotope(10*ones(2,1),0.5*diag(ones(2,1)));
U = zonotope([zeros(2,1),0.1*eye(2)]);

simOpts.x0 = randPoint(R0);
simOpts.u = randPoint(U,10);
simOpts.tFinal = 1;

[tCont,xCont] = simulate(linSys,simOpts);

% equivalent hybrid automaton
inv = polytope([1,0],1);
guard = conHyperplane([1,0],1);
reset.A = eye(2); reset.c = zeros(2,1);
trans = transition(guard,reset,2);
loc1 = location(inv,trans,linSys);

inv = polytope([-1,0],1);
guard = conHyperplane([1,0],-1);
reset.A = eye(2); reset.c = zeros(2,1);
trans = transition(guard,reset,1);
loc2 = location(inv,trans,linSys);

HA = hybridAutomaton([loc1;loc2]);

simOpts.startLoc = 2;

[tHyb,xHyb] = simulate(HA,simOpts);

tHyb = vertcat(tHyb{:});
xHyb = vertcat(xHyb{:});

% check if the two simulated trajectories are equivalent
xCont = interp1(tCont,xCont,tHyb);

if max(max(abs(xHyb-xCont))) > 0.1
    res(end+1,1) = false;
end


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
