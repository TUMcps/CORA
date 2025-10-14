function res = test_polygon_splitIntoConvexSets()
% test_polygon_splitIntoConvexSets - unit test function for 
%    polygon/splitIntoConvexSets
%
% Syntax:
%    res = test_polygon_splitIntoConvexSets()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Authors:       Tobias Ladner
% Written:       09-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% triangle
V_tria = [
    0 0 1;
    0 1 1
];
pgon = polygon(V_tria);
splits = pgon.splitIntoConvexSets();
assert(pgon == splits{1})

% stop sign
V_stop = [
    0 1 2 3 3 2 1 0;
    2 3 3 2 1 0 0 1 ;
] / 3;
pgon = polygon(V_stop);
splits = pgon.splitIntoConvexSets();
assert(pgon == splits{1})

% letter 'T'
V_T = [
    0 3 3 2 2 1 1 0
    3 3 2 2 0 0 2 2
];
pgon = polygon(V_T);
splits = pgon.splitIntoConvexSets();
assert(numel(splits) == 2)
assert(pgon == or(splits{1}, splits{2}))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
