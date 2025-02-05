function res = test_nn_rl_ctrlEnv_step()
% test_nn_rl_ctrlEnv_step - unit test function for
%     @ctrlEnvironment/ctrlEnvironment: Check given options are validated
%     correctly.
%
% Syntax:
%    res = test_nn_rl_ctrlEnv_step()
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
sys = nonlinearSys(f);
options.rl.env.x0 = interval([-4;0],[4;0]);
options.rl.env.initialOps = 'inf';
options.rl.env.dt = 0.1;
options.rl.timeStep = 0.02;

% Test ODE45
options.rl.env.solver = 'ODE45';
env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[env, observation] = env.reset();

action = single(1);
ops.u = action;
[~, x] = ode45(getfcn(env.ctrlDynamics, ops), single(0:options.rl.timeStep:options.rl.env.dt), [observation;action]);

[env, observation, ~, ~, ~] = env.step(action);

assert(all(x(end,1:2)' == observation) && env.stepNum == 2)

% Test Euler step
options.rl.env.solver = 'Euler';
env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);
[env, observation] = env.reset();

action = single(1);
ops.u = action;
t = options.rl.env.dt;
func = getfcn(env.ctrlDynamics, ops);
x = ([observation;action] + options.rl.env.dt*func(t,[observation;action]))';

[env, observation, ~, ~, ~] = env.step(action);

assert(all(x(end,1:2)' == observation) && env.stepNum == 2)

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Quadrocopter1D(x)
% reward function
if isnumeric(x)
    reward = -(abs(x(end,1))+0.01*abs(x(end,2)));
else
    w = [1,0.01,0];
    reward = -max(abs(supportFunc(x.timePoint.set{end},w,'lower')),abs(supportFunc(x.timePoint.set{end},w)));
end
end

function collisionBool = aux_collisionCheck_Quadrocopter1D(x)
% collision check
if isa(x,'numeric')
    if all(abs(x(:,1:2))<0.05,"all")
        collisionBool = true;
    else
        collisionBool = false;
    end
elseif isa(x,'reachSet')
    if norm(x.timePoint.set{1}.c)<0.05
        collisionBool = true;
    else
        collisionBool = false;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
