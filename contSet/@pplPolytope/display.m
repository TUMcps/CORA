function display(obj)
% display - Displays the C matrix and d vector of a pplPolytope
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
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display dimension
disp('dimension: ');
disp(obj.dim);

%display center
disp('C: ');
disp(obj.C);

%display generators
disp('d: ');
disp(obj.d); 

%------------- END OF CODE --------------