function display(obj)
% display - Displays a nonlinearSysDT object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - nonlinearSysDT object
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      27-October-2011
% Last update:  29-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

%display type
disp('type: Nonlinear time discrete system');

%display state space equations
%generate symbolic states
for i=1:obj.dim
    command=['syms x',num2str(i)];
    eval(command); 
    command=['x(',num2str(i),')=x',num2str(i),';'];
    eval(command);
end
%generate symbolic inputs
for i=1:obj.nrOfInputs
    command=['syms u',num2str(i)];
    eval(command); 
    command=['u(',num2str(i),')=u',num2str(i),';'];
    eval(command);
end

dx=obj.mFile(x,u)

disp('-----------------------------------');

%------------- END OF CODE --------------