function matP = matPolytope(intMat)
% matPolytope - converts an interval matrix to a matrix polytope
%
% Syntax:
%    matP = matPolytope(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
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

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   06-May-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%obtain vertices
V = vertices(intMat);

%instantiate matrix polytope
matP=matPolytope(V);

% ------------------------------ END OF CODE ------------------------------
