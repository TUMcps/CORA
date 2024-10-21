function display(nlnsys)
% display - displays a nonlinearSys object on the command window
%
% Syntax:
%    display(nlnsys)
%
% Inputs:
%    nlnsys - nonlinearSys object
%
% Outputs:
%    ---
%
% Example:
%    f = @(x,u) [x(2); ...
%               (1-x(1)^2)*x(2)-x(1)];
%    nlnsys = nonlinearSys('vanDerPol',f)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       17-October-2007
% Last update:   19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(nlnsys);

%display type
disp('Type: Nonlinear continuous-time system');

%create symbolic variables
vars = symVariables(nlnsys);

%insert symbolic variables into the system equations
f = nlnsys.mFile(vars.x,vars.u);

%display state space equations
disp('State-space equations:')
for i=1:length(f)
    disp(['  f(',num2str(i),') = ',char(f(i))]);
end

fprintf(newline);

%insert symbolic variables into the system equations
y = nlnsys.out_mFile(vars.x,vars.u);

% display output equation
disp('Output equations:')
for i=1:length(y)
    disp(['  y(',num2str(i),') = ',char(y(i))]);
end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
