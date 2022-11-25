function display(obj)
% display - Displays a nonlinDASys object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - nonlinDASys object
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
% Written:      27-October-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

% display number of algebric states
disp(['number of alg. variables: ', num2str(obj.nrOfConstraints)]);

%display type
disp('type: Nonlinear differential algebraic system');

%display state space equations
%create symbolic variables
var = symVariables(obj);

%insert symbolic variables into the system equations
fdyn=obj.dynFile(var.x,var.y,var.u);
disp('dynamic state space equations:')
for i=1:length(fdyn)
    disp(['f(',num2str(i),') = ',char(fdyn(i))]);
end

fcon=obj.conFile(var.x,var.y,var.u);
disp('constraint equations:')
for i=1:length(fcon)
    disp(['constr ',num2str(i),': 0 = ',char(fcon(i))]);
end

disp('-----------------------------------');


%------------- END OF CODE --------------