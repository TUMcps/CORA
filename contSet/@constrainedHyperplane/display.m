function display(obj)
% display - Displays the hyperplane and the constraint system
%
% Syntax:  
%    display(h)
%
% Inputs:
%    h - halfspace object
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
% Written:      10-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display hyperplane
display(obj.h);

%display constraint system
disp('Constraint sytem:')

disp('C:');
disp(obj.C);
disp('d:');
disp(obj.d);

%------------- END OF CODE --------------