function res = test_testCase
% test_testCase - unit test function for (deprecated) testCase 
%   constructor returning a trajectory object
%
% Syntax:
%    res = test_testCase()
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
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize 
nrStates = 4; nrOutputs = 2; nrInputs = 2; 
n_s = 10;
n_k = 8;
u = randn(n_k, nrInputs, n_s);
x = randn(n_k, nrStates, n_s);
y = randn(n_k, nrOutputs, n_s);
x0 = permute(x(1,:,:),[2 1 3]);
dt = 0.1;
name = "test";

% correct instantiations according to constructor
traj = testCase(y,u,x,dt);
assert(isa(traj, 'trajectory'))
traj = testCase(y,u(:,:,1),x,dt);
assert(isa(traj, 'trajectory'))
traj = testCase(y,u,x0,dt);
assert(isa(traj, 'trajectory'))
traj = testCase(y,u(:,:,1),x0(:,:,1),dt);
assert(isa(traj, 'trajectory'))
traj = testCase(y,u,x0,dt,name);
assert(isa(traj, 'trajectory'))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
