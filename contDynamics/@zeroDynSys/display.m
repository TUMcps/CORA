function display(obj)
% display - Displays a linearSys object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - linearSys object
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
% Written: 16-May-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

%display type
disp('type: Zero dynamics system');

disp('-----------------------------------');

%------------- END OF CODE --------------