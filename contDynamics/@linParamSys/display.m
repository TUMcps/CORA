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
% Written:      05-Aug-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

disp('-----------------------------------');

%display parent object
display@contDynamics(obj);

%display type
disp('type: Linear parameter varying system');

%display A-Matrix
disp('System matrix: ');
display(obj.A);

%display B-Matrix
disp('Input matrix: ');
display(obj.B);

%display Taylor terms
disp('Number of Taylor terms: ');
disp(obj.taylorTerms);

disp('-----------------------------------');

%------------- END OF CODE --------------