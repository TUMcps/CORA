function display(obj)
% display - Displays a location object
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - location object
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
% Written: 06-November-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display name
disp(['name: ',obj.name]);
    
%display continuous dynamics
display(obj.contDynamics);

disp('-----------------------------------');

%------------- END OF CODE --------------