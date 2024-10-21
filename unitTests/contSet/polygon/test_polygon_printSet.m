function res = test_polygon_printSet
% test_polygon_printSet - unit test function of printSet
%
% Syntax:
%    res = test_polygon_printSet
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
% See also: none

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test empty
pgon = polygon.empty(2);

printSet(pgon)
printSet(pgon,'high')
printSet(pgon,'high',true)
printSet(pgon,'high',false)

% test normal set
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
pgon = polygon(x(ind),y(ind));

printSet(pgon)
printSet(pgon,'high')
printSet(pgon,'high',true)
printSet(pgon,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
