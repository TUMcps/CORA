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

%display type
disp('type: Nonlinear differential algebraic system');

%display number of states, inputs and parameters
disp(['nr of states: ',num2str(obj.dim)]);
disp(['nr of inputs: ',num2str(obj.nrOfInputs)]);

%display state space equations
%create symbolic variables
[x,u]=symVariables(obj);

%insert symbolic variables into the system equations
t=0;
fdyn=obj.dynFile(t,x,u);
disp('dynamic state space equations:')
for i=1:length(fdyn)
    disp(['f(',num2str(i),') = ',char(fdyn(i))]);
end

fcon=obj.conFile(t,x,u);
disp('constraint equations:')
for i=1:length(fcon)
    disp(['constr ',num2str(i),': 0 = ',char(fcon(i))]);
end

disp('-----------------------------------');


%------------- END OF CODE --------------