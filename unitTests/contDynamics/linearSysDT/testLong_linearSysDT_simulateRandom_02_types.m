function res = testLong_linearSysDT_simulateRandom_02_types
% testLong_linearSysDT_simulateRandom_02_types - unit test function for simulating
%   different system dynamics
%
% Syntax:
%    res = testLong_linearSysDT_simulateRandom_02_types()
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

% linear system
A = [-0.3780    0.2839    0.5403   -0.2962
    0.1362    0.2742    0.5195    0.8266
    0.0502   -0.1051   -0.6572    0.3874
    1.0227   -0.4877    0.8342   -0.2372];
B = 0.25 * [-2 0 3;
    2 1 0;
    0 0 1;
    0 -2 1];
dt = 0.1;
sys = linearSysDT(linearSys(A,B),dt);
R0 = zonotope(randn(sys.nrOfDims, sys.nrOfDims+2));
U = zonotope(randn(sys.nrOfInputs, sys.nrOfInputs+2));

% simulate system with empty options
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
for type = ["standard", "gaussian", "rrt"] 
    options1 = options;
    options1.type = char(type);
    if type == "rrt"
        % compute reachable set
        optionsReach.timeStep = dt;
        optionsReach.taylorTerms = 4;
        optionsReach.zonotopeOrder = 20;
        R = reach(sys,params,optionsReach);

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
            paramsSim.u = traj(i).u;
            paramsSim.tFinal = traj(i).t(end);
            [t,x,~,y] = simulate(sys,paramsSim);
            assert(all(isapprox(x,traj(i).x,"tight"),'all'))
            assert(all(isapprox(t,traj(i).t,"tight"),'all'))
            assert(all(isapprox(y,traj(i).y,"tight"),'all'))
        end
    end
    % simulate
    traj = simulateRandom(sys,params,options1);

    % check correctness
    for i = 1:length(traj)
        paramsSim.x0 = traj(i).x(:,1);
        paramsSim.u = traj(i).u;
        paramsSim.tFinal = traj(i).t(end);
        [t,x,~,y] = simulate(sys,paramsSim);
        assert(all(isapprox(x,traj(i).x,"tight"),'all'))
        assert(all(isapprox(t,traj(i).t,"tight"),'all'))
        assert(all(isapprox(y,traj(i).y,"tight"),'all'))
    end
end

% gather results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
