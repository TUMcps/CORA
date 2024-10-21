function pgon = projectHighDim_(pgon,N,proj)
% projectHighDim_ - project a polygon to a higher-dimensional space
%
% Syntax:
%    pgon = projectHighDim_(pgon,N,proj)
%
% Inputs:
%    pgon - pgon object
%    N - dimension of the higher dimensional space
%    proj - states of the high dimensional space that correspond to the
%          states of the low dimensional polytope object
%
% Outputs:
%    pgon - polygon object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim, polygon/project

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if N ~= 2
    throw(CORAerror('CORA:wrongValue','second','Polygons are always two-dimensional.'))
end

% use project as polygons are always 2-dimensional
pgon = project(pgon,proj);

end

% ------------------------------ END OF CODE ------------------------------
