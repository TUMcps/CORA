function res = testLong_nonlinearSysDT_simulateRandom_01_types
% testLong_nonlinearSysDT_simulateRandom_01_types - unit test function for 
%   simulateRandom using different algorithms
%
% Syntax:
%    res = testLong_nonlinearSysDT_simulateRandom_01_types()
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
% Written:       10-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng(1)
n_k = 6;

% nonlinear system
dt = 0.01;
fun = @(x,u) x + dt*DOTBicycleDynamics_SRX_velEq(x,u);
sys = nonlinearSysDT('bicycle',fun, dt);
R0 = zonotope(randn(sys.nrOfDims, sys.nrOfDims+2));
U = zonotope(randn(sys.nrOfInputs, sys.nrOfInputs+2));

% simulate system with empty options
params.tStart = 0;
params.tFinal = (n_k-1)*dt;
params.U = U;
params.R0 = R0;
options = struct();
traj = simulateRandom(sys,params,options);

% define number of simulations
options.points = 3;
traj = simulateRandom(sys,params,options);
assert(length(traj) == 3)

% simulate different types
for type = ["standard", "gaussian","rrt"] 
    options1 = options;
    options1.type = char(type);
    dummy_u = [];
    if type == "rrt"
        dummy_u = zeros(sys.nrOfInputs,1); % required for simulate
        % use dummy reachable sets 
        R = struct();
        for i=1:n_k
            R.timePoint.time{i,1} = (i-1)*dt;
            R.timePoint.set{i,1} = zonotope(randn(sys.nrOfDims, sys.nrOfDims+2));
        end

        % set options
        options1.stretchFac = 1.2;
        options1.R = R;
        options1.vertSamp = true;

        % simulate
        traj = simulateRandom(sys,params,options1);
        options1.vertSamp = false;

        % check correctness
        for i = 1:length(traj)
            paramsSim.x0 = traj(i).x(:,1);
            paramsSim.u = [traj(i).u dummy_u];
            paramsSim.tFinal = traj(i).t(end);

            % simulate
            [t,x,~,y] = simulate(sys,paramsSim);
            assert(all(isapprox(x,traj(i).x,"tight"),'all'))
            assert(all(isapprox(t,traj(i).t,"tight"),'all'))
            %assert(all(isapprox(y,traj(i).y,"tight"),'all'))
        end
    end
    % simulate
    traj = simulateRandom(sys,params,options1);

    % check correctness
    for i = 1:length(traj)
        paramsSim.x0 = traj(i).x(:,1);
        paramsSim.u = [traj(i).u dummy_u];
        paramsSim.tFinal = traj(i).t(end);
        [t,x,~,y] = simulate(sys,paramsSim);
        assert(all(isapprox(x,traj(i).x,"tight"),'all'))
        assert(all(isapprox(t,traj(i).t,"tight"),'all'))
        %assert(all(isapprox(y,traj(i).y,"tight"),'all'))
    end
end

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
