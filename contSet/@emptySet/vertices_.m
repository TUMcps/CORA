function V = vertices_(O)
% vertices_ - returns the vertices of a empty set
%
% Syntax:
%    V = vertices_(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    V - vertices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/vertices

% Authors:       Tobier Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

V = vertices(dim(O),0);

% ------------------------------ END OF CODE ------------------------------
