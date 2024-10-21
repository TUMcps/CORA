function res = test_nonlinearSysDT_nonlinearSysDT
% test_nonlinearSysDT_nonlinearSysDT - unit test for constructor
%
% Syntax:
%    res = test_nonlinearSysDT_nonlinearSysDT
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

% Authors:       Mark Wetzlinger
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name
sysname = 'sys';
% sampling time
dt = 0.5;


% one-dimensional, without inputs
f_1D = @(x,u) x(1)^2;
sys = nonlinearSysDT(f_1D,dt);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_1D,dt);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_1D,dt,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(f_1D,dt,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)


% one-dimensional, with inputs
f_1D = @(x,u) x(1)^2 - u(1);
sys = nonlinearSysDT(f_1D,dt);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1 )
assert(sys.dt == dt)

sys = nonlinearSysDT(f_1D,dt,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_1D,dt,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

% three-dimensional
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-x(1); x(3)*x(2)];
sys = nonlinearSysDT(f_3D,dt);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 3)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_3D,dt,3,1);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 3)
assert(sys.dt == dt)

% with output equation
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-exp(x(1)); x(3)*x(2)];
g_1D = @(x,u) x(1)*x(2);
sys = nonlinearSysDT(f_3D,dt,g_1D);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_3D,dt,g_1D);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)

sys = nonlinearSysDT(sysname,f_3D,dt,3,1,g_1D,1);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)
assert(sys.dt == dt)


% wrong instantiations (should not reach next line)

% only name
assertThrowsAs(@nonlinearSysDT,'CORA:numInputArgsConstructor',sysname);

% name is numeric
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',1,f_1D,dt,1,1);
% only states, not inputs given
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',f_1D,dt,1);
% more states in output equation
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',f_1D,dt,g_1D);
% states and inputs, but not outputs given
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',f_3D,dt,3,1,g_1D);
% only x as input argument
f_x = @(x) x(1)^2;
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',f_x,dt);
% too many input arguments
f_xup = @(x,u,p) x(1)^2;
assertThrowsAs(@nonlinearSysDT,'CORA:wrongInputInConstructor',f_xup,dt);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
