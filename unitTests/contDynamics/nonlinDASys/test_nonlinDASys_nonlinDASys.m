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
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(sysname,f_1D,g_1D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 1 ...
        || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(f_1D,g_1D,1,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(sysname,f_1D,g_1D,1,1,1);
if ~strcmp(sys.name,sysname) || sys.dim ~= 1 ...
        || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end

% one-dimensional, with inputs
f_1D = @(x,y,u) x(1)^2 + u(1);
g_1D = @(x,y,u) x(1) - y(1);
sys = nonlinDASys(f_1D,g_1D);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(sysname,f_1D,g_1D);
if ~strcmp(sysname,sys.name) || sys.dim ~= 1 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(f_1D,g_1D,1,1,1);
if sys.dim ~= 1 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end
sys = nonlinDASys(sysname,f_1D,g_1D,1,1,1);
if ~strcmp(sysname,sys.name) || sys.dim ~= 1 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 1 || sys.nrOfOutputs ~= 1
    res = false;
end

% three-dimensional, no inputs, no constaints
f_3D = @(x,y,u) [sqrt(x(1)) - x(2)*u(1); x(2)-y(1); x(3)*x(2)];
g_2D = @(x,y,u) [x(1) - y(1); y(2) + x(2)];
sys = nonlinDASys(f_3D,g_2D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinDASys(f_3D,g_2D,3,1,2);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinDASys(sysname,f_3D,g_2D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end
sys = nonlinDASys(sysname,f_3D,g_2D,3,1,2);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 3
    res = false;
end

% with output equation
h_2D = @(x,y,u) [x(1)*y(1); x(2) - u(1)];
sys = nonlinDASys(f_3D,g_2D,h_2D);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinDASys(sysname,f_3D,g_2D,h_2D);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinDASys(f_3D,g_2D,3,1,2,h_2D,2);
if sys.dim ~= 3 || sys.nrOfInputs ~= 1 || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end
sys = nonlinDASys(sysname,f_3D,g_2D,3,1,2,h_2D,2);
if ~strcmp(sys.name,sysname) || sys.dim ~= 3 || sys.nrOfInputs ~= 1 ...
        || sys.nrOfConstraints ~= 2 || sys.nrOfOutputs ~= 2
    res = false;
end


% wrong instantiations

% only name
try
    nonlinDASys(sysname);
    res = false;
end
% name is numeric
try
    nonlinDASys(1,f_1D,g_1D,1,1);
    res = false;
end
% only states, not inputs or constraints given
try
    nonlinDASys(f_1D,g_1D,1);
    res = false;
end
% only states and inputs, but no constraints given
try
    nonlinDASys(f_1D,g_1D,1,1);
    res = false;
end
% more states in output equation
try
    nonlinDASys(f_1D,g_1D,h_2D);
    res = false;
end
% states, inputs, and constraints, but not outputs given
try
    nonlinDASys(f_3D,g_2D,3,1,2,h_2D);
    res = false;
end
% only dynamic function as input argument
try
    nonlinDASys(f_1D);
    res = false;
end
% only x as input argument
f_x = @(x) x(1)^2;
g_x = @(x) x(1);
try
    nonlinDASys(f_x,g_x);
    res = false;
end
% only x and y as input arguments
f_xy = @(x,y) x(1)^2 - y(1);
g_xy = @(x,y) x(1) + y(1);
try
    nonlinDASys(f_xy,g_xy);
    res = false;
end
% too many input arguments
f_xyup = @(x,y,u,p) x(1)^2;
try
    nonlinDASys(f_xyup);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
