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
%                25-April-2024 (TL, harmonized display with contSet classes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(matP)
    dispEmptyObj(matP,inputname(1));
    return
end


% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension, number of vertices
disp('dimension: ');
disp(dim(matP));
disp('nr of vertices: ');
disp(matP.numverts);

%display vertices
disp('vertices: ');
disp(matP.V)

% ------------------------------ END OF CODE ------------------------------
