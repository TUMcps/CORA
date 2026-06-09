function res = testLong_hybridAutomaton_falsifySpecial
% testLong_hybridAutomaton_falsifySpecial - unit test for falsification
%    featuring some special cases (nonlin. dynamics, nonlin. resets, etc.)
%
% Syntax:
%    res = testLong_hybridAutomaton_falsifySpecial
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

% Authors:       Niklas Kochdumper
% Written:       04-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% algorithm settings to test
alg = {'singleOpt','multiOpt'};
dyn = {'lin','mixInt'};


% Test 1: Nonlinear Dynamcis ----------------------------------------------

% load hybrid automaton
HA = poly3modes();

% modify the automaton by adding additive uncertainty as input
locs = [];

for i = 1:length(HA.location)
    loc = HA.location(i); sys = loc.contDynamics;
    sys = nonlinearSys(@(x,u) sys.mFile(x,u) + [u;0]);
    locs = [locs;location(loc.invariant,loc.transition,sys)];
end

HA = hybridAutomaton(locs);

% define reachability problem
params.R0 = interval([1;1],[1.5;1.5]);
params.U = zonotope(0.2*interval(-1,1)) + 0.05;
params.tFinal = 1.5;
params.startLoc = 3;

% % debug
% options_.alg = 'lin';
% options_.tensorOrder = 2;
% options_.timeStep = 0.01;
% options_.taylorTerms = 5;
% options_.zonotopeOrder = 20;
% options_.guardIntersect = 'polytope';
% options_.enclose = {'box'}; 
% options_.intersectInvariant = true;
% 
% R = reach(HA,params,options_);

% specification: unsafe set (falsifiable)
P = polytope([0 -1],-8.7);
spec = specification(P,'unsafeSet');
   
% loop over different algorithms
for i = 1:length(alg)
    
    options.falsifyAlg = alg{i};

    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == true);
    assert(any(contains(P,fals.traj.x)));
end


% Test 2: Multiple Modes --------------------------------------------------

% load hybrid automaton
HA = oscillator();

% modify the automaton by adding additive uncertainty as input
locs = [];

for i = 1:length(HA.location)

    loc = HA.location(i); sys = loc.contDynamics;
    sys = linearSys(sys.A,eye(2),sys.c);

    trans = [];

    for j = 1:length(loc.transition)
        trans_ = loc.transition(j);
        reset = linearReset(trans_.reset.A,eye(2),trans_.reset.c);
        trans = [trans;transition(trans_.guard,reset,trans_.target)];
    end

    locs = [locs;location(loc.invariant,trans,sys)];
end

HA = hybridAutomaton(locs);

% define reachability problem
params.R0 = zonotope([1;-1] + 0.1*interval(-ones(2,1),ones(2,1)));
params.U = zonotope(0.1*interval([-1;-1],[1;1])) + [0.05;-0.07];
params.tFinal = 3;
params.startLoc = 2;

% % debug
% options_.timeStep = 0.005;
% options_.taylorTerms = 5;
% options_.zonotopeOrder = 20;
% options_.guardIntersect = 'zonoGirard';
% options_.enclose = {'box','pca'}; 
% options_.intersectInvariant = true;
% 
% R = reach(HA,params,options_);

% specification: unsafe set (falsifiable)
P = polytope([-1 0],-0.55);
spec = specification(P,'unsafeSet',interval(2,3));
    
% loop over different algorithms
for i = 1:length(alg)
    for j = 1:length(dyn)
    
        options.falsifyAlg = alg{j};
        options.dynamics = dyn{j};
    
        [res,fals] = falsify(HA,params,options,spec);
        
        if strcmp(options.dynamics,'mixInt')
            assert(res == true);
            ind = find(contains(spec.time,fals.traj.t));
            assert(any(contains(P,fals.traj.x(:,ind))));
        end
    end
end


% Test 3: Nonlinear Resets ------------------------------------------------

% define hybrid automaton
f = @(x,u) [x(3); x(4); u(1); -9.81+0.01*x(4)^2];
sys = nonlinearSys(f);

syms x y vx vy
eq = -y + sin(x);
inv = levelSet(eq,[x;y;vx;vy],'<=');

guard = levelSet(-eq,[x;y;vx;vy],'==');

reset = nonlinearReset(...
    @(x,u) [x(1); ...
            x(2); ... 
            ((1-0.8*cos(x(1))^2)*x(3)+1.8*cos(x(1))*x(4))/(1+cos(x(1))^2); ...
            (1.8*cos(x(1))*x(3)+(-0.8+cos(x(1))^2)*x(4))/(1+cos(x(1))^2)]);

trans = transition(guard,reset,1);

loc = location(inv,trans,sys);

HA = hybridAutomaton(loc);

% define reachability problem
R0 = interval([0.48;4.98;0;-5],[0.52;5.02;0;-5]);

params.R0 = zonotope(R0); 
params.U = interval(-0.05,0.15);
params.startLoc = 1;                               
params.tFinal = 2.2;  

% % debug
% options_.alg = 'lin';
% options_.tensorOrder = 2;
% options_.timeStep = 0.01;
% options_.taylorTerms = 10;
% options_.zonotopeOrder = 20;
% options_.guardIntersect = 'levelSet';
% 
% R = reach(HA,params,options_);

% specification: unsafe set (falsifiable)
P = polytope([0 -1 0 0],-3.28);
spec = specification(P,'unsafeSet',interval(1,2));
    
% loop over different algorithms
for i = 1:length(alg)

    options.falsifyAlg = alg{i};

    [res,fals] = falsify(HA,params,spec);
        
    assert(res == true);
    ind = find(contains(spec.time,fals.traj.t));
    assert(any(contains(P,fals.traj.x(:,ind))));
end

% test successfull if it runs through without errors
res = true;

% ------------------------------ END OF CODE ------------------------------
