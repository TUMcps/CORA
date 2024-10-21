function pgon = enclosePoints(points)
% enclosePoints - enclose point cloud with polygon
%
% Syntax:
%    pgon = polygon.enclosePoints(points)
%
% Inputs:
%    points - numeric, vertices of polygon
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

ind = boundary(points', 0.5);
pgon = polygon(points(:, ind));

end

% ------------------------------ END OF CODE ------------------------------
