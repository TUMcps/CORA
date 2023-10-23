function display(sys)
% display - Displays a nonlinearSysDT object on the command window
%
% Syntax:
%    display(sys)
%
% Inputs:
%    sys - nonlinearSysDT object
%
% Outputs:
%    -
%
% Example:
%    f = @(x,u) [x(1) + u(1); x(2) + u(2)*cos(x(1)); x(3) + u(2)*sin(x(1))];
%    dt = 0.25;
%    sys = nonlinearSysDT(f,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       27-October-2011
% Last update:   29-January-2018 (NK)
%                19-June-2022 (MW)
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(sys);

%display type
disp('Type: Nonlinear discrete-time system');

% display sampling time
disp("Sampling time: " + sys.dt);

%create symbolic variables
vars = symVariables(sys);

%insert symbolic variables into the system equations
f = sys.mFile(vars.x,vars.u);

%display state space equations
disp('State-space equations:')
for i=1:length(f)
    disp(['  f(',num2str(i),') = ',char(f(i))]);
end

fprintf(newline);

%insert symbolic variables into the system equations
y = sys.out_mFile(vars.x,vars.u);

%display state space equations
disp('Output equations:')
for i=1:length(y)
    disp(['  y(',num2str(i),') = ',char(y(i))]);
end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
