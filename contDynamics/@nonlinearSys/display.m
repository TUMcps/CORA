function display(obj)
% display - Displays a nonlinearSys object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - nonlinearSys object
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

% Author: Matthias Althoff
% Written: 17-October-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

%display type
disp('type: Nonlinear time continuous system');

%display state space equations
x = sym('x',[obj.dim,1]);
u = sym('u',[obj.nrOfInputs,1]);

dx = obj.mFile(x,u)

disp('-----------------------------------');

%------------- END OF CODE --------------