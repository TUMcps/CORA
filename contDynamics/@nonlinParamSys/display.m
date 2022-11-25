function display(obj)
% display - Displays a linVarSys object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - linVarSys object
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-May-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

% display number of parameter
disp(['number of parameters: ',num2str(obj.nrOfParam)]);

%display type
disp('type: Nonlinear parameter system');

%display state space equations
%create symbolic variables
vars = symVariables(obj);

%insert symbolic variables into the system equations
f=obj.mFile(vars.x,vars.u,vars.p);
disp('state space equations:')
for i=1:length(f)
    disp(['f(',num2str(i),') = ',char(f(i))]);
end

disp('-----------------------------------');

%------------- END OF CODE --------------