function res = testLong_linearSys_reach_06_5dim_linAlg_all
% testLong_linearSys_reach_06_5dim_linAlg_all - unit_test_function of
%    linear reachability analysis with uncertain inputs (toy example),
%    all linear reach algorithms used (except krylov due to system size)
%
% Syntax:
%    res = testLong_linearSys_reach_06_5dim_linAlg_all
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
dim_x = length(A);
fiveDimSys = linearSys('fiveDimSys',A,B);

% model parameters
params.tFinal = 1;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));

% reachability settings
options.timeStep = 0.05;
options.taylorTerms = 4;
options.zonotopeOrder = 20;

% specification
spec = specification(polytope([0 -1 0 0 0],-2));

% reachability analysis
algs = {'standard','wrapping-free','fromStart','decomp'};
% note: krylov does not make sense for this small system dimension

Rs = cell(length(algs),1);
for i=1:length(algs)
	options.linAlg = algs{i};
    if strcmp(options.linAlg,'decomp')
        % additional options when decomp called
        options.partition = [1, 2; 3, 4; 5, 5];
        Rs{i} = reach(fiveDimSys, params, options, spec);
    else % 'standard', 'wrapping-free', 'fromStart'
        Rs{i} = reach(fiveDimSys, params, options, spec);
    end
end

% simulation
simOpt.points = 5;
traj = simulateRandom(fiveDimSys, params, simOpt);

% verify that end points of simulation are contain in reachable sets
finalPoints = cell2mat(arrayfun(@(s) s.x(:,end),traj,'UniformOutput',false)');
for i=1:numel(Rs)
    finalSet = Rs{i}.timePoint.set{end};
    assertLoop(all(contains(finalSet,finalPoints,'exact',1e-8)), i);
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
