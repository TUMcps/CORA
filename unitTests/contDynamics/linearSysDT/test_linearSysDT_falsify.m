function res = test_linearSysDT_falsify
% test_linearSysDT_falsify - unit test for falsification
%
% Syntax:
%    res = test_linearSysDT_falsify
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
% Written:       11-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % define linear discrete-time system
    A = [0.72 0.36; -0.18 1.08];
    dt = 0.2;

    sys = linearSysDT(A,eye(2),dt);
 
    % define reachability problem
    params.tFinal = 5;
    params.R0 = zonotope(interval([-10.1;9.9],[-9.9;10.1]));
    params.U = zonotope([-0.1;-0.1],[0.1;0.1]);

%     % debug
%     options_.zonotopeOrder = 50;
% 
%     R = reach(sys,params,options_);

    % loop over different algorithms
    alg = {'singleOpt','multiOpt','monteCarlo','koopman'};

    for i = 1:length(alg)

        options.falsifyAlg = alg{i};

        % Test 0: black-box algorithms
        if ismember(options.falsifyAlg,{'monteCarlo','koopman'})
            options.maxTime = 1;

            x = stl('x',2);
            phi = not(finally(globally(x(1) > 15,interval(0,0.4)),interval(0,3)));
            spec = specification(phi,'logic');
    
            P = interval([0;-2.5],[2;0]);
            spec = add(spec,specification(P,'unsafeSet'));
    
            falsify(sys,params,options,spec);

            continue;
        end

        % Test 1: unsafe set (falsifiable)
        P = interval([0;-2.5],[2;0]);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    
        % Test 2: unsafe set (not falsifiable)
        P = interval([-5;4.3],[0;6.3]);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        assert(all(~contains(P,fals.traj.x)));
    
        % Test 3: safe set (falsifiable)
        P = interval([-5;-5],[-2;-3]);
        spec = specification(P,'safeSet',interval(4.2));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 4: safe set (not falsifiable)
        P = interval([-5.1;-5],[-2;-3]);
        spec = specification(P,'safeSet',interval(4.2));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 5: multiple sets (falsifiable)
        P1 = interval([-5.1;-5],[-2;-3]);
        spec = specification(P1,'safeSet',interval(4.2));
    
        P2 = interval([0;-2.5],[2;0]);
        spec = add(spec,specification(P2,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P2,fals.traj.x)));
    
        % Test 6: multiple sets (not falsifiable)
        P1 = interval([-5.1;-5],[-2;-3]);
        spec = specification(P1,'safeSet',interval(4.2));
    
        P2 = interval([-5;4.3],[0;6.3]);
        spec = add(spec,specification(P2,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        [~,ind] = min((fals.traj.t - center(spec(1).time)).^2);
        assert(contains(P1,fals.traj.x(:,ind)));
        assert(~any(contains(P2,fals.traj.x)));
    
        % Test 7: temporal logic (falsifiable)
        x = stl('x',2);
        phi = not(finally(globally(x(1) > 15,interval(0,0.1)),interval(0,3)));
    
        [res,fals] = falsify(sys,params,options,phi);
    
        assert(res == true);
        ind = find(fals.traj.x(1,:) > 15);
        assert(max(fals.traj.t(ind)) - min(fals.traj.t(ind)) > 0.1);
    
        % Test 8: temporal logic (not falsifiable)
        x = stl('x',2);
        phi = not(finally(globally(x(1) > 15,interval(0,0.4)),interval(0,3)));
    
        [res,fals] = falsify(sys,params,options,phi);
    
        assert(res == false);
        ind = find(fals.traj.x(1,:) > 15);
        assert(max(fals.traj.t(ind)) - min(fals.traj.t(ind)) < 0.4);
    
        % Test 9: system with output equation
        sys_ = linearSys(sys.A,sys.B,sys.c,[1,1]);
    
        P = polytope(-1,-27);
        spec = specification(P,'unsafeSet',interval(0,2));
    
        [res,fals] = falsify(sys_,params,options,spec);
    
        assert(res == true);
        ind = find(fals.traj.t <= 2);
        assert(any(fals.traj.y(ind) > 27));
    
        % Test 10: system with time-varying input
        params_ = params;
        params_.u = [2 -2 3 -1 2; -1 -2 2 1 0];
    
        P = polytope([-1,0],-11);
        spec = specification(P,'unsafeSet',interval(2,5));
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        ind = find(fals.traj.t > 2);
        assert(any(contains(P,fals.traj.x(:,ind))));
    
        % Test 11: non-zero initial time
        params_ = params;
        params_.tStart = 0.5;
        params_.tFinal = 5.5;
    
        P = interval([-5;-5],[-2;-3]);
        spec = specification(P,'safeSet',interval(4.2) + 0.5);
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));

        % Test 12: mixed reach-avoid and temporal logic specification
        x = stl('x',2);
        phi = not(finally(globally(x(1) > 15,interval(0,0.4)),interval(0,3)));
        spec = specification(phi,'logic');

        P = interval([0;-2.5],[2;0]);
        spec = add(spec,specification(P,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));

        % Test 13: system with constant input
        sys_ = linearSysDT(sys.A,sys.B,[1;-2],sys.dt);

        P = polytope([1,0],-27.8);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys_,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    end

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
