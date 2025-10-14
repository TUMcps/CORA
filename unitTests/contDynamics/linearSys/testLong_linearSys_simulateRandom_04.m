function res = testLong_linearSys_simulateRandom_04
% testLong_linearSys_simulateRandom_04 - unit test for simulation of linearSys,
%    with different simulation types
%
% Syntax:
%    res = testLong_linearSys_simulateRandom_04
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

% Authors:       Tobias Ladner
% Written:       22-May-2024
% Last update:   10-September-2025 (LL, check correctness of results)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Model Parameters --------------------------------------------------------

tFinal = 5;
params = {};
% case: only R0
params.tFinal = tFinal;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);
params.U = zonotope(interval([0.9; -0.25; -0.1], [1.1; 0.25; 0.1]));

% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
dim_x = length(A);          % n
dim_u = dim(params.U);   % m
B = randn(dim_x,dim_u);
sys = linearSys('fiveDimSys',A,B);

% Random Simulations ------------------------------------------------------

N = 10;
simOpt.points = N;

% Check for completion ----------------------------------------------------

% compute reachable set
optionsReach.timeStep = 0.2;
optionsReach.taylorTerms = 4;
optionsReach.zonotopeOrder = 20;
R = reach(sys,params,optionsReach);

types = {'standard','gaussian','rrt'}; %

for j=1:numel(types)
    % generate seed
    seed = randi(2^32);
    rng(seed);
    simOpt.type = types{j};

    % set options
    simOpt_ = simOpt;
    if strcmp(types{j}, 'rrt')
        simOpt_.stretchFac = 1.2;
        simOpt_.R = R;
        simOpt_.vertSamp = true;
    end
    traj = simulateRandom(sys, params, simOpt_);
    assert(all(size(traj) == [N,1]),"seed: %s",int2str(seed));

    % check correctness
    for i = 1:length(traj)
        paramsSim.x0 = traj(i).x(:,1);
        paramsSim.u = traj(i).u;
        paramsSim.tFinal = traj(i).t(end);
        [t,x,~,y] = simulate(sys,paramsSim);
        assert(all(isapprox(x(:,1),traj(i).x(:,1),"tight"),'all'))
        assert(all(isapprox(t(:,1),traj(i).t(:,1),"tight"),'all'))
        %assert(all(isapprox(x(:,end),traj(i).x(:,end),"tight"),'all'))
        assert(all(isapprox(t(:,end),traj(i).t(:,end),"tight"),'all'))
    end
end

% combine results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
