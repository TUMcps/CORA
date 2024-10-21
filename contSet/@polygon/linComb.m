function pgon = linComb(pgon1, pgon2)
% linComb - compute linear combination of two polygons
%
% Syntax:
%    pgon = linComb(pgon1, pgon2)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


w = warning();
warning('off');

list1 = triangulation(pgon1);
list2 = triangulation(pgon2);

pgon = [];

for i = 1:length(list1)
    for j = 1:length(list2)
        pgon = pgon | convHull(list1{i}, list2{j});
    end
end

warning(w);
end

% ------------------------------ END OF CODE ------------------------------
