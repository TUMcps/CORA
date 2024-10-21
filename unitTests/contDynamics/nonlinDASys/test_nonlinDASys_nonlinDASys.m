function res = test_nonlinDASys_nonlinDASys
% test_nonlinDASys_nonlinDASys - unit test for constructor
%
% Syntax:
%    res = test_nonlinDASys_nonlinDASys
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

% assume true, wait for failure
res = true;

% name
sysname = 'sys';


% one-dimensional, no inputs
f_1D = @(x,y,u) x(1)^2;
g_1D = @(x,y,u) y(1) + x(1);
sys = nonlinDASys(f_1D,g_1D);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(sysname,f_1D,g_1D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(f_1D,g_1D,1,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(sysname,f_1D,g_1D,1,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)


% one-dimensional, with inputs
f_1D = @(x,y,u) x(1)^2 + u(1);
g_1D = @(x,y,u) x(1) - y(1);
sys = nonlinDASys(f_1D,g_1D);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(sysname,f_1D,g_1D);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(f_1D,g_1D,1,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinDASys(sysname,f_1D,g_1D,1,1,1);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 1)
assert(sys.nrOfOutputs == 1)

% three-dimensional, no inputs, no constaints
f_3D = @(x,y,u) [sqrt(x(1)) - x(2)*u(1); x(2)-y(1); x(3)*x(2)];
g_2D = @(x,y,u) [x(1) - y(1); y(2) + x(2)];
sys = nonlinDASys(f_3D,g_2D);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinDASys(f_3D,g_2D,3,1,2);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinDASys(sysname,f_3D,g_2D);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinDASys(sysname,f_3D,g_2D,3,1,2);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 3)

% with output equation
h_2D = @(x,y,u) [x(1)*y(1); x(2) - u(1)];
sys = nonlinDASys(f_3D,g_2D,h_2D);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinDASys(sysname,f_3D,g_2D,h_2D);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinDASys(f_3D,g_2D,3,1,2,h_2D,2);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinDASys(sysname,f_3D,g_2D,3,1,2,h_2D,2);
assert(strcmp(sysname,sys.name))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfConstraints == 2)
assert(sys.nrOfOutputs == 2)


% wrong instantiations

% only name
assertThrowsAs(@nonlinDASys,'CORA:numInputArgsConstructor',sysname);
% name is numeric
assertThrowsAs(@nonlinDASys,'MATLAB:UndefinedFunction',1,f_1D,g_1D,1,1);
% only states, not inputs or constraints given
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_1D,g_1D,1);
% only states and inputs, but no constraints given
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_1D,g_1D,1,1);
% more states in output equation
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_1D,g_1D,h_2D);
% states, inputs, and constraints, but not outputs given
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_3D,g_2D,3,1,2,h_2D);
% only dynamic function as input argument
assertThrowsAs(@nonlinDASys,'CORA:numInputArgsConstructor',f_1D);
% only x as input argument
f_x = @(x) x(1)^2;
g_x = @(x) x(1);
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_x,g_x);
% only x and y as input arguments
f_xy = @(x,y) x(1)^2 - y(1);
g_xy = @(x,y) x(1) + y(1);
assertThrowsAs(@nonlinDASys,'CORA:wrongInputInConstructor',f_xy,g_xy);
% too many input arguments
f_xyup = @(x,y,u,p) x(1)^2;
assertThrowsAs(@nonlinDASys,'CORA:numInputArgsConstructor',f_xyup);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
