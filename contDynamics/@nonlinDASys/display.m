function display(nlnsysDA)
% display - displays a nonlinDASys object on the command window
%
% Syntax:
%    display(nlnsysDA)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%
% Outputs:
%    -
%
% Example:
%    f = @(x,y,u) x(1)+1+u(1);
%    g = @(x,y,u) (x(1)+1)*y(1) + 2;
%    nlnsysDA = nonlinDASys(f,g)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       27-October-2011
% Last update:   19-June-2022
%                23-November-2022 (TL, dispInput)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(nlnsysDA);

%display type
disp("Type: Nonlinear differential-algebraic system");

% display number of algebric states
disp("  number of alg. variables: " + nlnsysDA.nrOfConstraints);

%create symbolic variables
vars = symVariables(nlnsysDA);

fprintf(newline);

% insert symbolic variables into the dynamic function
fdyn=nlnsysDA.dynFile(vars.x,vars.y,vars.u);
disp('State-space equations: f(x,y,u)')
for i=1:length(fdyn)
    disp(['  f(',num2str(i),') = ',char(fdyn(i))]);
end

fprintf(newline);

% insert symbolic variables into the constraint function
fcon=nlnsysDA.conFile(vars.x,vars.y,vars.u);
disp('Constraint equations: g(x,y,u)')
for i=1:length(fcon)
    disp(['  constraint ',num2str(i),': 0 = ',char(fcon(i))]);
end

fprintf(newline);

% insert symbolic variables into the constraint function
fout=nlnsysDA.out_mFile(vars.x,vars.y,vars.u);
disp('Output equations: h(x,y,u)')
for i=1:length(fout)
    disp(['  h(',num2str(i),') = ',char(fout(i))]);
end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
