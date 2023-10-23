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

% assume true, wait for failure
res = true;

% name
sysname = 'sys';


% one-dimensional, no inputs, no parameters
f_1D = @(x,u,p) x(1)^2;
sys = nonlinParamSys(f_1D);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(sysname,f_1D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 1 ...
        || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(f_1D,1,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(sysname,f_1D,1,1,1);
if ~strcmp(sys.name,sysname) || sys.dim ~= 1 ...
        || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end

% one-dimensional, with inputs, with parameters
f_1D = @(x,u,p) x(1)^2 + u(1)*p(1);
sys = nonlinParamSys(f_1D);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(sysname,f_1D);
if ~strcmp(sysname,sys.name) || sys.dim ~= 1 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(f_1D,1,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinParamSys(sysname,f_1D,1,1,1);
if ~strcmp(sysname,sys.name) || sys.dim ~= 1 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end

% three-dimensional, no inputs, no parameters
f_3D = @(x,u,p) [sqrt(x(1)) - x(2)*u(1); x(2)-p(1); x(3)*p(2)];
sys = nonlinParamSys(f_3D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinParamSys(f_3D,3,1,2);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinParamSys(sysname,f_3D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinParamSys(sysname,f_3D,3,1,2);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end

% with output equation
g_2D = @(x,u,p) [x(1)*p(1); x(2) - u(1)];
sys = nonlinParamSys(f_3D,g_2D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinParamSys(sysname,f_3D,g_2D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinParamSys(f_3D,3,1,2,g_2D,2);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinParamSys(sysname,f_3D,3,1,2,g_2D,2);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfParam ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end


% wrong instantiations
if CHECKS_ENABLED

% only name
try
    nonlinParamSys(sysname);
    res = false;
end
% name is numeric
try
    nonlinParamSys(1,f_1D,1,1);
    res = false;
end
% only states, not inputs or params given
try
    nonlinParamSys(f_1D,1);
    res = false;
end
% only states and inputs, but no params given
try
    nonlinParamSys(f_1D,1,1);
    res = false;
end
% more states in output equation
try
    nonlinParamSys(f_1D,g_2D);
    res = false;
end
% states, inputs, and params, but not outputs given
try
    nonlinParamSys(f_3D,3,1,2,g_2D);
    res = false;
end
% only x as input argument
f_x = @(x) x(1)^2;
try
    nonlinParamSys(f_x);
    res = false;
end
% only x and u as input arguments
f_xu = @(x,u) x(1)^2 - u(1);
try
    nonlinParamSys(f_xu);
    res = false;
end
% too many input arguments
f_xyup = @(x,y,u,p) x(1)^2;
try
    nonlinParamSys(f_xyup);
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
