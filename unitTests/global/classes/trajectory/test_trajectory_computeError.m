function res = test_trajectory_computeError
% test_trajectory_computeError - unit test function for computing the
%   error between system simulations and the real data
%
% Syntax:
%    res = test_trajectory_computeError()
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

% Authors:       Laura Luetzow
% Written:       29-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng(1)
n_k = 6;

for dynamics = ["pedestrian", "lorenz", "linearSys", "nonlinearSys"]
    % define system
    if dynamics == "linearSys"
        % linear system
        A = [-0.3780    0.2839    0.5403   -0.2962
            0.1362    0.2742    0.5195    0.8266
            0.0502   -0.1051   -0.6572    0.3874
            1.0227   -0.4877    0.8342   -0.2372];
        B = 0.25 * [-2 0 3;
            2 1 0;
            0 0 1;
            0 -2 1];
        sys = linearSys(A,B);
        dt = 0.1;
    elseif dynamics == "nonlinearSys"
        % nonlinear system
        fun = @(x,u) [x(2); (1-x(1)^2)*x(2)-x(1)];
        sys = nonlinearSys('vanDerPol',fun);
        dt = 0.1;
    else
        sys = loadDynamics(dynamics);
        dt = sys.dt;
    end

    % sample initial state and inputs
    x0_1 = randn(sys.nrOfDims,1);
    u1 = randn(sys.nrOfInputs,n_k);
    x0_2 = randn(sys.nrOfDims,1);
    u2 = randn(sys.nrOfInputs,n_k);

    % simulate system
    params.x0 = x0_1;
    if isa(sys,'nonlinearSysDT') 
        params.u = u1(:,1:end-1);
    else
        params.u = u1;
    end
    params.tStart = 0;
    params.tFinal = (n_k-1)*dt;
    t = 0:dt:(n_k-1)*dt;
    [t1,x1] = simulate(sys,params);
    % interpolate states
    if isa(sys,'nonlinearSys') || isa(sys,'linearSys')
        [~,ind] = unique(t1);
        x1 = interp1(t1(ind),x1(:,ind)',t,'linear','extrap')';
        t1 = t;
    end

    % simulate again
    params.x0 = x0_2;
    if isa(sys,'nonlinearSysDT')
        params.u = u2(:,1:end-1);
    else
        params.u = u2;
    end
    [t2,x2] = simulate(sys,params);
    % interpolate states
    if isa(sys,'nonlinearSys') || isa(sys,'linearSys')
        [~,ind] = unique(t2);
        x2 = interp1(t2(ind),x2(:,ind)',t,'linear','extrap')';
        t2 = t;
    end

    % create trajectories
    traj1 = trajectory(u1,x1,[],t1);
    traj2 = trajectory(u2,x2,[],t2);
    traj = [traj1; traj2];

    % compute error (should be approximately zero)
    err12 = computeError(traj,sys);
    assert(abs(err12) < 1e-6);

    % create random trajectories
    x3 = randn(sys.nrOfDims,length(t2));
    u3 = randn(sys.nrOfInputs,length(t2));
    traj3 = trajectory(u3,x3,[],t2);

    % compute error (should be bigger than 0)
    err3 = computeError(traj3,sys);
    assert(abs(err3) > 1e-3);

    % compute error for all trajectories
    err123 = computeError([traj1; traj2; traj3],sys);
    assert(abs(err123-err12-err3) < 1e-6);
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
