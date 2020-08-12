function display(obj)
% display - Displays the C matrix and d vector of a mptPolytope
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - pplPolytope object
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

%display mpt polytope
display(obj.P); 
fprintf(newline);

%------------- END OF CODE --------------