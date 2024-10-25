function pgon = expandBoundaries(pgon, d)
% expandBoundaries - expands the boundaries of a polygon by a distance d
%
% Syntax:
%    pgon = expandBoundaries(pgon,d)
%
% Inputs:
%    pgon - polygon
%    d - numeric, distance to expand the boundaries
%
% Outputs:
%    pgon - expanded polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pgon.set = polybuffer(pgon.set, d);
pgon.V = [];

end

% ------------------------------ END OF CODE ------------------------------
