function pgon = project(pgon,dims)
% project - projects a polygon onto the specified dimensions
%
% Syntax:
%    pgon = project(pgon,dims)
%
% Inputs:
%    pgon - polygon
%    dims - dimensions for projection
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check inputs
inputArgsCheck({ ...
  {pgon,'att','polygon'}, ...
  {dims,'att','numeric','nonnegative'}, ...
});
if numel(dims) > 2 || max(dims) > 2
    throw(CORAerror('CORA:wrongValue','second','Polygons can only be two-dimensional.'))
end

% consider 1D
if isscalar(dims)
    % convert to interval
    pgon = project(interval(pgon),dims);
    return
end

% get vertices
V = vertices_(pgon);

% project (dims can only be [1 1],[1 2],[2 1],[2 2])
V = V(dims,:);

% obtain polygon result
pgon = polygon(V);

end

% ------------------------------ END OF CODE ------------------------------
