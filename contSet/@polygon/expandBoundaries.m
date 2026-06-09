function pgon = expandBoundaries(pgon, d, varargin)
% expandBoundaries - expands the boundaries of a polygon by a distance d
%
% Syntax:
%    pgon = expandBoundaries(pgon,d,varargin)
%
% Inputs:
%    pgon - polygon
%    d - numeric, distance to expand the boundaries
%    varargin - additional arguments for polyshape/polybuffer
%
% Outputs:
%    pgon - expanded polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyshape/polybuffer

% Authors:       Tobias Ladner
% Written:       08-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pgon.set = polybuffer(pgon.set, d, varargin{:});
pgon.V = [];

end

% ------------------------------ END OF CODE ------------------------------
