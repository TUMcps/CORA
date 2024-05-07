function V = vertices(matP)
% vertices - returns the vertices of a matrix polytope
%
% Syntax:
%    V = vertices(matP)
%
% Inputs:
%    matP - matPolytope object
%
% Outputs:
%    V - numeric (n x m x N)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       24-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%return vertices
V=matP.V;

% ------------------------------ END OF CODE ------------------------------
