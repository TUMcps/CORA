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

disp('-----------------------------------');

%------------- END OF CODE --------------