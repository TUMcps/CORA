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

% assume true, wait for failure
res = true;

% name
sysname = 'sys';
% sampling time
dt = 0.5;


% one-dimensional, without inputs
f_1D = @(x,u) x(1)^2;
sys = nonlinearSysDT(f_1D,dt);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_1D,dt);
if ~strcmp(sys.name,sysname) ...
        || sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_1D,dt,1,1);
if ~strcmp(sys.name,sysname) ...
        || sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(f_1D,dt,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end

% one-dimensional, with inputs
f_1D = @(x,u) x(1)^2 - u(1);
sys = nonlinearSysDT(f_1D,dt);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(f_1D,dt,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_1D,dt,1,1);
if ~strcmp(sys.name,sysname) || sys.dim ~= 1 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfOutputs ~= 1 || sys.dt ~= dt
    res = false;
end


% three-dimensional
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-x(1); x(3)*x(2)];
sys = nonlinearSysDT(f_3D,dt);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 3 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_3D,dt,3,1);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 3 ...
        || sys.dt ~= dt
    res = false;
end

% with output equation
f_3D = @(x,u) [sqrt(x(1)) - x(2)*u(1); x(2)-exp(x(1)); x(3)*x(2)];
g_1D = @(x,u) x(1)*x(2);
sys = nonlinearSysDT(f_3D,dt,g_1D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_3D,dt,g_1D);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end
sys = nonlinearSysDT(sysname,f_3D,dt,3,1,g_1D,1);
if ~strcmp(sysname,sys.name) ...
        || sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfOutputs ~= 1 ...
        || sys.dt ~= dt
    res = false;
end


% wrong instantiations (should not reach next line)

% only name
try
    nonlinearSysDT(sysname);
    res = false;
end
% name is numeric
try
    nonlinearSysDT(1,f_1D,dt,1,1);
    res = false;
end
% only states, not inputs given
try
    nonlinearSysDT(f_1D,dt,1);
    res = false;
end
% more states in output equation
try
    nonlinearSysDT(f_1D,dt,g_1D);
    res = false;
end
% states and inputs, but not outputs given
try
    nonlinearSysDT(f_3D,dt,3,1,g_1D);
    res = false;
end
% only x as input argument
f_x = @(x) x(1)^2;
try
    nonlinearSysDT(f_x,dt);
    res = false;
end
% too many input arguments
f_xup = @(x,u,p) x(1)^2;
try
    nonlinearSysDT(f_xup,dt);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
