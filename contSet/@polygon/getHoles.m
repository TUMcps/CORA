function holes = getHoles(pgon)
% getHoles - returns the holes of the polygon
%
% Syntax:
%    pgons = getHoles(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    holes - cell of polygons
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyshape/holes

% Authors:       Tobias Ladner
% Written:       13-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

holes = arrayfun(@(pshape) polygon(pshape), pgon.set.holes,'UniformOutput',false);

end

% ------------------------------ END OF CODE ------------------------------
