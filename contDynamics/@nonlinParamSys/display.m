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

%display type
disp('type: Nonlinear parameter system');

%display number of states, inputs and parameters
disp(['nr of states: ',num2str(obj.dim)]);
disp(['nr of inputs: ',num2str(obj.nrOfInputs)]);
disp(['nr of parameters: ',num2str(obj.nrOfParam)]);

%display state space equations
%create symbolic variables
[x,u,~,p]=symVariables(obj);

%insert symbolic variables into the system equations
t=0;
f=obj.mFile(t,x,u,p);
disp('state space equations:')
for i=1:length(f)
    disp(['f(',num2str(i),') = ',char(f(i))]);
end

disp('-----------------------------------');

%------------- END OF CODE --------------