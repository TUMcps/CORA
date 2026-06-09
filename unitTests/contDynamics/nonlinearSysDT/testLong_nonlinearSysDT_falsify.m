function res = testLong_nonlinearSysDT_falsify
% testLong_nonlinearSysDT_falsify - unit test for falsification
%
% Syntax:
%    res = testLong_nonlinearSysDT_falsify
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
% Written:       12-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % define nonlinear system
    dt = 0.5;
    fun = @(x,u) tank6EqDT(x,u,dt);

    sys = nonlinearSysDT(fun,dt);
 
    % define reachability problem
    params.tFinal = 400;
    params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);
    params.U = zonotope([0,0.005]);

%     % debug
%     options_.taylorTerms = 4;
%     options_.zonotopeOrder = 50;
%     options_.alg = 'lin';
%     options_.tensorOrder = 2;
%     options_.errorOrder = 10;
% 
%     options_.lagrangeRem.simplify = 'simplify';
% 
%     R = reach(sys,params,options_);

    % loop over different algorithms
    alg = {'singleOpt','multiOpt','monteCarlo','koopman'};

    for i = 1:length(alg)

        options.falsifyAlg = alg{i};
        options.maxTime = 60;

        % Test 0: black-box algorithms
        if ismember(options.falsifyAlg,{'monteCarlo','koopman'})
            options.maxTime = 1;

            P = polytope([0 0 0 1 0 0],1.78);
            x = stl('x',6);
            phi = not(finally(globally(x(1) < 1.6,interval(0,100)), ...
                                                        interval(0,300)));
            spec = specification(P,'unsafeSet');
            spec = add(spec,specification(phi,'logic'));
    
            falsify(sys,params,options,spec);

            continue;
        end

        % Test 1: unsafe set (falsifiable)
        P = polytope([0 0 0 1 0 0],1.78);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    
        % Test 2: unsafe set (not falsifiable)
        P = polytope([0 0 0 1 0 0],1.5);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        assert(all(~contains(P,fals.traj.x)));
    
        % Test 3: safe set (not falsifiable)
        I = interval([2.4;2.4;-10*ones(4,1)],[3.4;3.4;10*ones(4,1)]);
        spec = specification(I,'safeSet',interval(282));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == false);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 4: safe set (falsifiable)
        I = interval([2.3;2.12;-10*ones(4,1)],[3.2;2.8;10*ones(4,1)]);
        spec = specification(I,'safeSet',interval(282));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 5: multiple sets (falsifiable)
        I = interval([2.4;2.4;-10*ones(4,1)],[3.4;3.4;10*ones(4,1)]);
        spec = specification(I,'safeSet',interval(282));
        
        P = polytope([0 0 0 1 0 0],1.78);
        spec = add(spec,specification(P,'unsafeSet'));
        
        [res,fals] = falsify(sys,params,options,spec);
        
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
        
        % Test 6: multiple sets (not falsifiable)
        I = interval([2.4;2.4;-10*ones(4,1)],[3.4;3.4;10*ones(4,1)]);
        spec = specification(I,'safeSet',interval(282));
        
        P = polytope([0 0 0 1 0 0],1.5);
        spec = add(spec,specification(P,'unsafeSet'));
        
        [res,fals] = falsify(sys,params,options,spec);
        
        assert(res == false);
        assert(all(~contains(P,fals.traj.x)));
    
        % Test 7: temporal logic (not falsifiable)
        x = stl('x',6);
        phi = not(finally(globally(x(1) < 1.4,interval(0,50)),interval(0,300)));
        
        [res,fals] = falsify(sys,params,options,phi);
    
        assert(res == false);
        ind = find(fals.traj.x(1,:) < 1.4);
        assert(isempty(ind) || max(fals.traj.t(ind)) - min(fals.traj.t(ind)) < 50);
    
        % Test 8: temporal logic (falsifiable)
        x = stl('x',6);
        phi = not(finally(globally(x(1) < 1.6,interval(0,50)),interval(0,300)));
        
        [res,fals] = falsify(sys,params,options,phi);
    
        assert(res == true);
        ind = find(fals.traj.x(1,:) < 1.6);
        assert(max(fals.traj.t(ind)) - min(fals.traj.t(ind)) > 50);
    
        % Test 9: system with output equation
        sys_ = nonlinearSys(@tank6Eq,@(x,u) x(4)^2);
    
        P = polytope(1,1.8^2);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys_,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.y)));
    
        % Test 10: system with time-varying inputs
        params_ = params;
        params_.u = [-0.004 -0.002 -0.001 0.003];
    
        P = polytope([0 0 0 1 0 0],1.7);
        spec = specification(P,'unsafeSet');
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    
        % Test 11: non-zero initial time
        params_ = params;
        params_.tStart = 100;
        params_.tFinal = 500;
    
        I = interval([2.3;2.12;-10*ones(4,1)],[3.2;2.8;10*ones(4,1)]);
        spec = specification(I,'safeSet',interval(382));
    
        [res,fals] = falsify(sys,params_,options,spec);
    
        assert(res == true);
        [~,ind] = min((fals.traj.t - center(spec.time)).^2);
        assert(~contains(spec.set,fals.traj.x(:,ind)));
    
        % Test 12: mixed reach-avoid and temporal logic specification
        x = stl('x',6);
        phi = not(finally(globally(x(1) < 1.4,interval(0,50)),interval(0,300)));
        spec = specification(phi,'logic');
    
        P = polytope([0 0 0 1 0 0],1.78);
        spec = add(spec,specification(P,'unsafeSet'));
    
        [res,fals] = falsify(sys,params,options,spec);
    
        assert(res == true);
        assert(any(contains(P,fals.traj.x)));
    end

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
