function res = test_nn_rl_ctrlEnv_validateOptions()
% test_nn_rl_ctrlEnv_validateOptions - unit test function for
%     @ctrlEnvironment/ctrlEnvironment: Check given options are validated
%     correctly.
%
% Syntax:
%    res = test_nn_rl_ctrlEnv_validateOptions()
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
% See also:

% Authors:       Manuel Wendl
% Written:       27-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------
m = 0.05;
g = 9.81;
% open-loop system
f = @(x,u) [x(2);(u(1)+1)/(2*m)-g];
sys = contDynamics('Name',1,2);

% test wrong values
assertThrowsAs(@ctrlEnvironment,'CORA:wrongValue',...
    sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D);

sys = nonlinearSys(f);
options.rl.env.x0 = [[-4;0],[4;0]];

assertThrowsAs(@ctrlEnvironment,'CORA:wrongFieldValue',...
    sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);

% test ctrlEnvironment.reset
options.rl.env.x0 = interval([-4;0],[4;0]);
options.rl.env.initialOps = 'None';

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[~, observation] = env.reset();
% Check with no initial ops
assert(all(observation == [0;0]))

options.rl.env.initialOps = 'inf';

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[~, observation] = env.reset();

% Check with inf initial ops
assert(all(observation == options.rl.env.x0.inf))

options.rl.env.initialOps = 'sup';

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[~, observation] = env.reset();

% Check with sup initial ops
assert(all(observation == options.rl.env.x0.sup))

% test outOfDomain
options.rl.env.dt = .1;
options.rl.env.timeStep = .2;

assertThrowsAs(@ctrlEnvironment,'CORA:outOfDomain',...
    sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);

options.rl.env.timeStep = .02;
options.rl.env.maxSteps = 0;

assertThrowsAs(@ctrlEnvironment,'CORA:outOfDomain',...
    sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);

% test reward function
options.rl.env.maxSteps = 30;
options.rl.env.evalMode = 'point';

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[~, observation] = env.reset();

assert(isnumeric(observation))

options.rl.env.evalMode = 'set';

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[~, observation] = env.reset();

assert(isa(observation,'zonotope'))

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Quadrocopter1D(x)
% reward function of quadrocopter

if isnumeric(x)
    reward = -(abs(x(end,1))+0.01*abs(x(end,2)));
else
    w = [1,0.01,0];
    reward = -max(abs(supportFunc(x.timePoint.set{end},w,'lower')),...
        abs(supportFunc(x.timePoint.set{end},w)));
end

end

function collisionBool = aux_collisionCheck_Quadrocopter1D(x)
% collision check
if isa(x,'numeric')
    collisionBool = all(abs(x(:,1:2)) < 0.05,"all");
elseif isa(x,'reachSet')
    collisionBool = norm(x.timePoint.set{1}.c) < 0.05;
end

end

% ------------------------------ END OF CODE ------------------------------
