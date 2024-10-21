function res = test_nonlinParamSys_nonlinParamSys
% test_nonlinParamSys_nonlinParamSys - unit test for constructor
%
% Syntax:
%    res = test_nonlinParamSys_nonlinParamSys
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


% one-dimensional, no inputs, no parameters
f_1D = @(x,u,p) x(1)^2;
sys = nonlinParamSys(f_1D);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(sysname,f_1D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(f_1D,1,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(sysname,f_1D,1,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)


% one-dimensional, with inputs, with parameters
f_1D = @(x,u,p) x(1)^2 + u(1)*p(1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(sysname,f_1D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(f_1D,1,1,1);
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

sys = nonlinParamSys(sysname,f_1D,1,1,1);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 1)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 1)
assert(sys.nrOfOutputs == 1)

% three-dimensional, no inputs, no parameters
f_3D = @(x,u,p) [sqrt(x(1)) - x(2)*u(1); x(2)-p(1); x(3)*p(2)];
sys = nonlinParamSys(f_3D);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinParamSys(f_3D,3,1,2);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinParamSys(sysname,f_3D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 3)

sys = nonlinParamSys(sysname,f_3D,3,1,2);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 3)

% with output equation
g_2D = @(x,u,p) [x(1)*p(1); x(2) - u(1)];
sys = nonlinParamSys(f_3D,g_2D);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinParamSys(sysname,f_3D,g_2D);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinParamSys(f_3D,3,1,2,g_2D,2);
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 2)

sys = nonlinParamSys(sysname,f_3D,3,1,2,g_2D,2);
assert(strcmp(sys.name,sysname))
assert(sys.nrOfStates == 3)
assert(sys.nrOfInputs == 1)
assert(sys.nrOfParam == 2)
assert(sys.nrOfOutputs == 2)


% only name
assertThrowsAs(@nonlinParamSys,'MATLAB:UndefinedFunction',sysname);
% name is numeric
assertThrowsAs(@nonlinParamSys,'MATLAB:UndefinedFunction',1,f_1D,1,1);
% only states and inputs, but no params given
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_1D,1,1);
% more states in output equation
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_1D,g_2D);
% states, inputs, and params, but not outputs given
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_3D,3,1,2,g_2D);
% only x as input argument
f_x = @(x) x(1)^2;
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_x);
% only x and u as input arguments
f_xu = @(x,u) x(1)^2 - u(1);
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_xu);
% too many input arguments
f_xyup = @(x,y,u,p) x(1)^2;
assertThrowsAs(@nonlinParamSys,'CORA:wrongInputInConstructor',f_xyup);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
