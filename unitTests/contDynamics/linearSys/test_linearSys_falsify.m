function res = test_linearSys_falsify
% test_linearSys_falsify - unit test for falsification
%
% Syntax:
%    res = test_linearSys_falsify
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

    % define linear system
    sys = linearSys([-0.7 -2; 2 -0.7],1);
 
    % define reachability problem
    params.tFinal = 1;
    params.R0 = zonotope(interval([9.9;9.9],[10.1;10.1]));
    params.U = zonotope([-0.5;-0.5],[0.5;0.5]);

%     % debug
%     options.timeStep = 0.005;
%     options.taylorTerms = 10;
%     options.zonotopeOrder = 50;
% 
%     R = reach(sys,params,options);

    % loop over different algorithms
    alg = {'singleOpt','multiOpt','monteCarlo','koopman'};

    for i = 1:length(alg)

        options.falsifyAlg = alg{i};

        % Test 0: black-box algorithms
        if ismember(options.falsifyAlg,{'monteCarlo','koopman'})
            options.maxTime = 1;

            P = interval([-5.1;4.3],[0;6.3]);
            x = stl('x',2);
            phi = until(x(2) > 4,x(1) < -6.2,interval(0.5,1));
            spec = specification(P,'unsafeSet');
            spec = add(spec,specification(phi,'logic'));
    
            falsify(sys,params,options,spec);

            continue;
        end

        % Test 1: unsafe set (falsifiable)
        P = interval([-5.1;4.3],[0;6.3]);
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
        P = interval([-2.4;9],[-2;9.8]);
        spec = specification(P,'safeSet',interval(0.505));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 4: safe set (not falsifiable)
        P = interval([-2.6;8.9],[-2;10]);
        spec = specification(P,'safeSet',interval(0.505));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 5: multiple sets (falsifiable)
        P1 = interval([-2.6;8.9],[-2;10]);
        spec = specification(P1,'safeSet',interval(0.505));
    
        P2 = interval([-5.1;4.3],[0;6.3]);
        spec = add(spec,specification(P2,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P2,fals.traj.x)));
    
        % Test 6: multiple sets (not falsifiable)
        P1 = interval([-2.6;8.9],[-2;10]);
        spec = specification(P1,'safeSet',interval(0.505));
    
        P2 = interval([-5;4.3],[0;6.3]);
        spec = add(spec,specification(P2,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        [~,ind] = min((fals.traj.t - center(spec(1).time)).^2);
        assert(contains(P1,fals.traj.x(:,ind)));
        assert(~any(contains(P2,fals.traj.x)));
    
        % Test 7: temporal logic (falsifiable)
        x = stl('x',2);
        phi1 = finally(globally(x(2) > 11,interval(0,0.2)),interval(0,0.8));
        phi2 = until(x(2) > 4,x(1) < -6.2,interval(0.5,1));
    
        [res1,fals1] = falsify(sys,params,options,phi1);
        [res2,fals2] = falsify(sys,params,options,phi2);
    
        assert(res1 == true);
        ind = find(fals1.traj.x(2,:) > 11);
        assert(max(fals1.traj.t(ind)) - min(fals1.traj.t(ind)) < 0.2);
    
        assert(res2 == true);
        ind = find(fals2.traj.t >= 0.5 & fals2.traj.t <= 1);
        assert(any(contains(polytope([-1,0;0 1],[6.2;4]),fals2.traj.x(:,ind))));
    
        % Test 8: temporal logic (not falsifiable)
        x = stl('x',2);
        phi1 = finally(globally(x(2) > 11,interval(0,0.08)),interval(0,0.9));
        phi2 = not(until(x(1) > -6,x(2) < 3.9,interval(0.5,1)));
    
        [res1,fals1] = falsify(sys,params,options,phi1);
        [res2,fals2] = falsify(sys,params,options,phi2);
    
        assert(res1 == false);
        ind = find(fals1.traj.x(2,:) > 11);
        assert(max(fals1.traj.t(ind)) - min(fals1.traj.t(ind)) > 0.08);

        assert(res2 == false);
        ind = find(fals2.traj.t >= 0.5 & fals2.traj.t <= 1);
        assert(~any(contains(polytope([-1,0;0 1],[6;3.9]),fals2.traj.x(:,ind))));
    
        % Test 9: system with output equation
        sys_ = linearSys(sys.A,sys.B,sys.c,[1,1]);
    
        P = polytope(1,10);
        spec = specification(P,'unsafeSet',interval(0,0.4));
    
        [res,fals] = falsify(sys_,params,options,spec);
    
        assert(res == true);
        ind = find(fals.traj.t <= 0.4);
        assert(any(fals.traj.y(ind) < 10));
    
        % Test 10: system with time-varying input
        params_ = params;
        params_.u = [2 -2 3 -1; -1 -2 2 1];
    
        P = interval([-4.8;5.2],[4;6]);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    
        % Test 11: non-zero initial time
        params_ = params;
        params_.tStart = 0.5;
        params_.tFinal = 1.5;
    
        P = interval([-2.4;9],[-2;9.8]);
        spec = specification(P,'safeSet',interval(0.505) + 0.5);
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));

        % Test 12: mixed reach-avoid and temporal logic specification
        x = stl('x',2);
        phi = finally(globally(x(2) > 11,interval(0,0.08)),interval(0,0.9));
        spec = specification(phi,'logic');

        P = interval([-5.1;4.3],[0;6.3]);
        spec = add(spec,specification(P,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));

        % Test 13: system with constant input
        sys_ = linearSys(sys.A,sys.B,[1;-2]);

        P = interval([-4;5],[-3;6]);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys_,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    end

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
