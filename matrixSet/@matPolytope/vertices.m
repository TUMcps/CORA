function matV = vertices(matP)
% vertices - returns the vertices of a matrix polytope
%
% Syntax:  
%    matV = vertices(matP)
%
% Inputs:
%    matP - matrix polytope
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

% Author:       Matthias Althoff
% Written:      24-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%return vertices
matV=matP.vertex;


%------------- END OF CODE --------------