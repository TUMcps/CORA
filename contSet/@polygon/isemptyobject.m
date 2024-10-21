function res = isemptyobject(pgon)
% isemptyobject - checks if the polygon is entirely empty
%
% Syntax:
%    res = isemptyobject(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, use vertices_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isempty(vertices_(pgon));

end

% ------------------------------ END OF CODE ------------------------------
