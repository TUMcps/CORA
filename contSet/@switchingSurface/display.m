function display(obj)
% display - Displays the normal vector and distance to the origin of a
% halfspace
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
% Written:      06-June-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display dimension, generators
disp('dimension: ');
disp(length(obj.c));

%display normal vector
disp('normal vector: ');
disp(obj.c);

%display distance to origin
disp('distance to origin: ');
disp(obj.d);

%------------- END OF CODE --------------