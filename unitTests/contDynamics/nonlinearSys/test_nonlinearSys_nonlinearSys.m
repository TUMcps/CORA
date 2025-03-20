function res = test_nonlinearSys_nonlinearSys
% test_nonlinearSys_nonlinearSys - unit test for constructor
%
% Syntax:
%    res = test_nonlinearSys_nonlinearSys
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
% Last update:   10-November-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true, wait for failure
res = true;

% name
sysname = 'sys';

% one-dimensional, without inputs
f_1D = @(x,u) x(1)^2;
sys = nonlinearSys(f_1D);
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

% one-dimensional with name
sys = nonlinearSys(sysname,f_1D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

% one-dimensional, name, dims specified
sys = nonlinearSys(sysname,f_1D,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

% one-dimensional, dims specified
sys = nonlinearSys(f_1D,1,1);
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

% one-dimensional, with inputs
f_1D = @(x,u) x(1)^2 - u(1);
sys = nonlinearSys(f_1D);
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinearSys(f_1D,1,1);
% with dims specified
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinearSys(sysname,f_1D,1,1);
% with name/dims specified
assert(strcmp(sys.name,sysname))
assert(sys.nrOfDims == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfOutputs == 1)

% three-dimensional
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-x(1); x(3)*x(2)];
sys = nonlinearSys(f_3D);
% multi-dimensional, no further arguments
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 3);

sys = nonlinearSys(sysname,f_3D,3,1);
% multi-dimensional, dims specified
assert(strcmp(sysname,sys.name))
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 3);

% explicitly state that there are no inputs, but CORA needs at least one
% input internally...
sys = nonlinearSys(sysname,f_3D,3,0);
% multi-dimensional, dims, name specified
assert(strcmp(sysname,sys.name))
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 3);

% with output equation
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-exp(x(1)); x(3)*x(2)];
g_1D = @(x,u) x(1)*x(2);
sys = nonlinearSys(f_3D,g_1D);
% multi-dimensional, output and no other parameters
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 1);

sys = nonlinearSys(sysname,f_3D,g_1D);
% multi-dimensional, output and name specified
assert(strcmp(sysname,sys.name))
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 1);

sys = nonlinearSys(sysname,f_3D,3,1,g_1D,1);
% multi-dimensional, output, dims and name specified
assert(strcmp(sysname,sys.name))
assert(sys.nrOfDims == 3);
assert(sys.nrOfInputs == 1);
assert(sys.nrOfOutputs == 1);


% wrong instantiations

% only name
assertThrowsAs(@nonlinearSys,'MATLAB:UndefinedFunction',sysname);

% name is numeric
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',1,f_1D,1,1);

% only states, not inputs given
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',f_1D,1);

% more states in output equation
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',f_1D,g_1D);

% states and inputs, but not outputs given
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',f_3D,3,1,g_1D);

% only x as input argument
f_x = @(x) x(1)^2;
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',f_x);

% too many input arguments
f_xup = @(x,u,p) x(1)^2;
assertThrowsAs(@nonlinearSys,'CORA:wrongInputInConstructor',f_xup);

% ------------------------------ END OF CODE ------------------------------
