function matP = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or a 
%    matrix polytope with a matrix polytope
%
% Syntax:
%    matP = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or matrix polytope
%    factor2 - numerical matrix or matrix polytope
%
% Outputs:
%    matP - matrix polytope
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix polytope
    matP=factor2;
    %compute vertices
    for i=1:matP.verts
        matP.vertex{i}=matrix*matP.vertex{i};
    end
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix polytope
    matP=factor1;
    %compute vertices
    for i=1:matP.verts
        matP.vertex{i}=matP.vertex{i}*matrix;
    end
    
%both factors are polytope matrices
else
    %initialize matrix zonotope
    matP1=factor1;
    %initialize matrix zonotope
    matP2=factor2;
    %initialize matrix zonotope
    matP=matPolytope();
    %compute vertices
    for j=1:matP1.verts
        for i=1:matP2.verts
            matP.vertex{end+1}=matP1.vertex{j}*matP2.vertex{i};
        end
    end
    %update number of generators
    matP.verts=matP1.verts*matP2.verts;
end

% ------------------------------ END OF CODE ------------------------------
