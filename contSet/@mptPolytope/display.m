function display(P)
% display - Displays the properties of a mptPolytope object (constraint
%    matrix and offset vector) on the command window
%
% Syntax:  
%    display(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  09-June-2020 (MW, add variable name)
% Last revision:---

%------------- BEGIN CODE --------------

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display polytope
display(P.P); 
fprintf(newline);

%------------- END OF CODE --------------