function matP = matPolytope(matI)
% matPolytopes - converts an interval matrix to a matrix polytope
%
% Syntax:  
%    matP = matPolytope(matI)
%
% Inputs:
%    matI - interval matrix 
%
% Outputs:
%    matP - polytope matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain vertices
V = dominantVertices(matI,inf);

%instantiate matrix polytope
matP=matPolytope(V);

%------------- END OF CODE --------------