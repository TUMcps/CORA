function res = test_linearSysDT_identify
% test_linearSysDT_identify - unit test for system identification
%
% Syntax:
%    res = test_linearSysDT_identify
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
% See also: linearSysDT/identify

% Authors:       Niklas Kochdumper, Laura Luetzow
% Written:       17-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % define algorithms to test
    alg = {'dmd','opt'};

    % loop over all algorithms
    for i = 1:length(alg)
        for j = 1:2
    
            % generate random test data
            dt = rand(); n = randi([2,4]); m = randi([1,3]); 
            N1 = randi([3,10]); N2 = randi([3,10]); s = randi([1,10]);
    
            t1 = (0:N1-1)*dt; t2 = (0:N2-1)*dt; t = {t1;t2};
            x1 = s*(rand(n,N1)-0.5); x2 = s*(rand(n,N2)-0.5); x = {x1;x2};
            u1 = s*(rand(m,N1)-0.5); u2 = s*(rand(m,N2)-0.5); u = {u1;u2};
    
            % set algorithm options
            options = [];
            options.alg = alg{i};

            if j == 2
                s = 0.45+rand()*0.55;
                options.dt = s*dt;
            end

            % Test 1: single trajectory without input
            sys = linearSysDT.identify(x1,t1,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 2: single trajectory with input
            sys = linearSysDT.identify(x1,t1,u1,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 3: multiple trajectories without input
            sys = linearSysDT.identify(x,t,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 4: multiple trajectories with input
            sys = linearSysDT.identify(x,t,u,options);
            assert(isa(sys,'linearSysDT'));

            % Test 5: single trajectory object without input
            traj = trajectory([],x1,[],t1);
            sys = linearSysDT.identify(traj,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 6: single trajectory object with input
            traj = trajectory(u1,x1,[],t1); 
            sys = linearSysDT.identify(traj,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 7: multiple trajectory objects without input
            traj = [trajectory([],x1,[],t1); trajectory([],x2,[],t2)];
            sys = linearSysDT.identify(traj,options);
            assert(isa(sys,'linearSysDT'));
    
            % Test 8: multiple trajectory objects with input
            traj = [trajectory(u1,x1,[],t1); trajectory(u2,x2,[],t2)];
            sys = linearSysDT.identify(traj,options);
            assert(isa(sys,'linearSysDT'));
        end
    end

    % test successfull if it runs through without errors
    res = true;

% ------------------------------ END OF CODE ------------------------------
