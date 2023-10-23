function display(matP)
% display - Displays the vertices of a matrix polytope
%
% Syntax:
%    display(matP)
%
% Inputs:
%    matP - matPolytope object
%
% Outputs:
%    -
%
% Example:
%    matP = matPolytope()
%    matP = matPolytope({[2 3; 2 1],[3 4; 3 2],[1 1; 1 0]})
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-June-2010
% Last update:   03-April-2023 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(matP)
    
    dispEmptyObj(matP,inputname(1));

else

    %display dimension, number of vertices
    disp('dimension: ');
    disp(dim(matP));
    disp('nr of vertices: ');
    disp(matP.verts);
    
    %display vertices
    disp('vertices: ');
    for i=1:length(matP.vertex)
        disp(matP.vertex{i}); 
        disp('---------------'); 
    end
end

% ------------------------------ END OF CODE ------------------------------
