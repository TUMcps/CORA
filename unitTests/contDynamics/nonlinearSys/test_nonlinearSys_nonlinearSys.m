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
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(sysname,f_1D);
if ~strcmp(sys.name,sysname) ...
        || sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(sysname,f_1D,1,1);
if ~strcmp(sys.name,sysname) ...
        || sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(f_1D,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end

% one-dimensional, with inputs
f_1D = @(x,u) x(1)^2 - u(1);
sys = nonlinearSys(f_1D);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(f_1D,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(sysname,f_1D,1,1);
if ~strcmp(sys.name,sysname) ...
        || sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end


% three-dimensional
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-x(1); x(3)*x(2)];
sys = nonlinearSys(f_3D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinearSys(sysname,f_3D,3,1);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 3
    res = false;
end
% explicitly state that there are no inputs, but CORA needs at least one
% input internally...
sys = nonlinearSys(sysname,f_3D,3,0);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 3
    res = false;
end

% with output equation
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-exp(x(1)); x(3)*x(2)];
g_1D = @(x,u) x(1)*x(2);
sys = nonlinearSys(f_3D,g_1D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(sysname,f_3D,g_1D);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinearSys(sysname,f_3D,3,1,g_1D,1);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end


% wrong instantiations (should not reach next line)
if CHECKS_ENABLED

% only name
try
    nonlinearSys(sysname);
    res = false;
end
% name is numeric
try
    nonlinearSys(1,f_1D,1,1);
    res = false;
end
% only states, not inputs given
try
    nonlinearSys(f_1D,1);
    res = false;
end
% more states in output equation
try
    nonlinearSys(f_1D,g_1D);
    res = false;
end
% states and inputs, but not outputs given
try
    nonlinearSys(f_3D,3,1,g_1D);
    res = false;
end
% only x as input argument
f_x = @(x) x(1)^2;
try
    nonlinearSys(f_x);
    res = false;
end
% too many input arguments
f_xup = @(x,u,p) x(1)^2;
try
    nonlinearSys(f_xup);
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
