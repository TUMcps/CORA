function P = polytope(O)
% polytope - converts a emptySet to a polytope object
%
% Syntax:
%    P = polytope(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    O = emptySet(2);
%    P = polytope(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/empty

% Authors:       Tobias Ladner
% Written:       25-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

P = polytope.empty(O.dimension);

% ------------------------------ END OF CODE ------------------------------
