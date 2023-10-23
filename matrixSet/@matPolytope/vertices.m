function matV = vertices(matP)
% vertices - returns the vertices of a matrix polytope
%
% Syntax:
%    matV = vertices(matP)
%
% Inputs:
%    matP - matPolytope object
%
% Outputs:
%    matV - cell array of matrix vertices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       24-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%return vertices
matV=matP.vertex;


% ------------------------------ END OF CODE ------------------------------
