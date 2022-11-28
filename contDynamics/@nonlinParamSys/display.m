function display(sys)
% display - Displays a nonlinParamSys object on the command window
%
% Syntax:  
%    display(sys)
%
% Inputs:
%    sys - nonlinParamSys object
%
% Outputs:
%    ---
%
% Example:
%    f = @(x,u,p) [x(2); ...
%                  p(1)*(1-x(1)^2)*x(2)-x(1)];
%    g = @(x,u,p) x(1);
%    sys = nonlinParamSys('vanDerPol',f,g)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      27-May-2011
% Last update:  19-June-2022
%               23-November-2022 (TL: dispInput)
% Last revision:---

%------------- BEGIN CODE --------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(sys);

%display type
disp("Type: Nonlinear continuous-time parametric system");

% display number of parameters
disp("  number of parameters: " + sys.nrOfParam);
fprintf(newline);

%create symbolic variables
vars = symVariables(sys);

%insert symbolic variables into the system equations
f = sys.mFile(vars.x,vars.u,vars.p);

%display state space equations
disp('State-space equations:')
for i=1:length(f)
    disp(['  f(',num2str(i),') = ',char(f(i))]);
end

fprintf(newline);

%insert symbolic variables into the system equations
y = sys.out_mFile(vars.x,vars.u,vars.p);

%display state space equations
disp('Output equations:')
for i=1:length(y)
    disp(['  y(',num2str(i),') = ',char(y(i))]);
end

fprintf(newline);

%------------- END OF CODE --------------