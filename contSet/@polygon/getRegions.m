function regions = getRegions(pgon)
% getRegions - returns the regions of the polygon
%
% Syntax:
%    pgons = getRegions(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    regions - cell of polygons
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyshape/regions

% Authors:       Tobias Ladner
% Written:       13-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

regions = arrayfun(@(pshape) polygon(pshape), pgon.set.regions,'UniformOutput',false);

end

% ------------------------------ END OF CODE ------------------------------
