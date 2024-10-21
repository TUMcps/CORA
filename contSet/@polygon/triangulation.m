function list = triangulation(pgon)
% triangulation - compute an triangulation of the polygon
%
% Syntax:
%    list = triangulation(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    list - cell of polygon triangles
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute triangulation
T = triangulation(pgon.set);

% convert the triangles to polygons
list = cell(size(T.ConnectivityList, 1), 1);

for i = 1:size(T.ConnectivityList, 1)
    V = T.Points(T.ConnectivityList(i, :), :);
    list{i} = polygon(V(:, 1), V(:, 2));
end

end

% ------------------------------ END OF CODE ------------------------------
