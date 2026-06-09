function res = testLong_hybridAutomaton_falsify
% testLong_hybridAutomaton_falsify - unit test for falsification
%
% Syntax:
%    res = testLong_hybridAutomaton_falsify
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
% Written:       27-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load hybrid automaton
HA = bouncing_ball(-0.9);

% modify the automaton by adding additive uncertainty as input
locs = [];

for i = 1:length(HA.location)
    loc = HA.location(i); sys = loc.contDynamics;
    sys = linearSys(sys.A,[0;1],sys.c);
    locs = [locs;location(loc.invariant,loc.transition,sys)];
end

HA = hybridAutomaton(locs);

% define reachability problem
params.R0 = zonotope([1;0] + 0.1*interval(-ones(2,1),ones(2,1)));
params.U = zonotope(0.1*interval(-1,1)) + 0.05;
params.tFinal = 1;
params.startLoc = 1;

% % debug
% options_.timeStep = 0.005;
% options_.taylorTerms = 10;
% options_.zonotopeOrder = 40;
% options_.guardIntersect = 'polytope';
% options_.enclose = {'box'}; 
% options_.intersectInvariant = true;
% 
% R = reach(HA,params,options_);

% loop over different algorithms
alg = {'singleOpt','multiOpt','monteCarlo','koopman'};

for i = 1:length(alg)

    options.falsifyAlg = alg{i};

    % Test 0: black-box algorithms
    if ismember(options.falsifyAlg,{'monteCarlo','koopman'})
        options.maxTime = 1;

        x = stl('x',2);
        phi = not(finally(globally(x(1) > 0.9,interval(0,0.3)),interval(0,0.7)));
        spec = specification(phi,'logic');
    
        P = polytope([-1 0],-0.9);
        spec = add(spec,specification(P,'unsafeSet',interval(0.4,1)));
    
        falsify(HA,params,options,spec);

        continue;
    end

    % Test 1: unsafe set (falsifiable)
    P = polytope([-1 0],-0.9);
    spec = specification(P,'unsafeSet',interval(0.4,1));
        
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == true);
    ind = find(contains(spec.time,fals.traj.t));
    assert(any(contains(P,fals.traj.x(:,ind))));
    
    % Test 2: unsafe set (not falsifiable)
    P = polytope([-1 0],-0.98);
    spec = specification(P,'unsafeSet',interval(0.4,1));
        
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == false);
    ind = find(contains(spec.time,fals.traj.t));
    assert(~any(contains(P,fals.traj.x(:,ind))));
    
    % Test 3: safe set (not falsifiable)
    I = interval([0.37;1.5],[0.7;3.2]);
    spec = specification(I,'safeSet',interval(0.62));
    
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == false);
    [~,ind] = min(abs(fals.traj.t - center(spec.time)));
    assert(contains(I,fals.traj.x(:,ind)));
    
    % Test 4: safe set (falsifiable)
    I = interval([0.37;1.5],[0.57;3.2]);
    spec = specification(I,'safeSet',interval(0.62));
    
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == true);
    [~,ind] = min(abs(fals.traj.t - center(spec.time)));
    assert(~contains(I,fals.traj.x(:,ind)));
    
    % Test 5: multiple sets (falsifiable)
    P = polytope([-1 0],-0.98);
    spec = specification(P,'unsafeSet',interval(0.4,1));
    
    I = interval([0.37;1.5],[0.54;3.2]);
    spec = add(spec,specification(I,'safeSet',interval(0.62)));
    
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == true);
    [~,ind] = min(abs(fals.traj.t - center(spec(2).time)));
    assert(~contains(I,fals.traj.x(:,ind)));
    
    % Test 6: multiple sets (not falsifiable)
    P = polytope([-1 0],-0.98);
    spec = specification(P,'unsafeSet',interval(0.4,1));
    
    I = interval([0.37;1.5],[0.7;3.2]);
    spec = add(spec,specification(I,'safeSet',interval(0.62)));
    
    [res,fals] = falsify(HA,params,options,spec);
        
    assert(res == false);
    ind = find(contains(spec(1).time,fals.traj.t));
    assert(~any(contains(P,fals.traj.x(:,ind))));
    [~,ind] = min(abs(fals.traj.t - center(spec(2).time)));
    assert(contains(I,fals.traj.x(:,ind)));
    
    % Test 7: temporal logic (not falsifiable)
    x = stl('x',2);
    phi = not(finally(globally(x(1) > 0.9,interval(0,0.3)),interval(0,0.7)));
    
    [res,fals] = falsify(HA,params,options,phi);
        
    assert(res == false);
    ind = find(fals.traj.x(1,:) > 0.9  & fals.traj.t > 0.5);
    assert(isempty(ind) || max(fals.traj.t(ind)) - min(fals.traj.t(ind)) < 0.3);
    
    % Test 8: temporal logic (falsifiable)
    x = stl('x',2);
    phi = not(finally(globally(x(1) > 0.9,interval(0,0.05)),interval(0.5,0.9)));
    
    [res,fals] = falsify(HA,params,options,phi);
        
    assert(res == true);
    ind = find(fals.traj.x(1,:) > 0.9  & fals.traj.t > 0.5);
    assert(max(fals.traj.t(ind)) - min(fals.traj.t(ind)) > 0.05);
    
    % Test 9: system with time-varying inputs
    params_ = params;
    params_.u = [-0.1 -0.1 0.08 0.08];
    
    P = polytope([-1 0],-0.91);
    spec = specification(P,'unsafeSet',interval(0.4,1));
        
    [res,fals] = falsify(HA,params_,options,spec);
        
    assert(res == true);
    ind = find(contains(spec.time,fals.traj.t));
    assert(any(contains(P,fals.traj.x(:,ind))));
    
    % Test 10: non-zero initial time
    params_ = params;
    params_.tStart = 0.5;
    params_.tFinal = 1.5;
    
    I = interval([0.37;1.5],[0.57;3.2]);
    spec = specification(I,'safeSet',interval(1.12));
    
    [res,fals] = falsify(HA,params_,options,spec);
        
    assert(res == true);
    [~,ind] = min(abs(fals.traj.t - center(spec.time)));
    assert(~contains(I,fals.traj.x(:,ind)));
        
    % Test 11: mixed reach-avoid and temporal logic specification
    x = stl('x',2);
    phi = not(finally(globally(x(1) > 0.9,interval(0,0.3)),interval(0,0.7)));
    spec = specification(phi,'logic');
    
    P = polytope([-1 0],-0.9);
    spec = add(spec,specification(P,'unsafeSet',interval(0.4,1)));
    
    [res,fals] = falsify(HA,params,options,spec);
    
    assert(res == true);
    ind = find(contains(spec(2).time,fals.traj.t));
    assert(any(contains(P,fals.traj.x(:,ind))));
end

% test successfull if it runs through without errors
res = true;


% ------------------------------ END OF CODE ------------------------------
