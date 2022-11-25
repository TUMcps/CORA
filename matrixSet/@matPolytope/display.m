function display(matP)
% display - Displays the vertices of a matrix polytope
%
% Syntax:  
%    display(matP)
%
% Inputs:
%    matP - matrix polytope
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
% Written:      21-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display dimension, generators
disp('dimension: ');
disp(matP.dim);
disp('nr of vertices: ');
disp(matP.verts);

%display vertices
disp('vertices: ');
for i=1:length(matP.vertex)
    disp(matP.vertex{i}); 
    disp('---------------'); 
end

%------------- END OF CODE --------------