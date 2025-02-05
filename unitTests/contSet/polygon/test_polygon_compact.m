function res = test_polygon_compact()
% test_polygon_compact - unit test function for polygon/compact
%
% Syntax:
%    res = test_polygon_compact()
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

% Authors:       Niklas Kochdumper
% Written:       02-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);

% test 'douglasPeuker'
pgon_ = compact(pgon,'douglasPeucker',0.1);

assert(contains(pgon_,pgon));

% test 'simplify'
pgon_ = compact(pgon,'simplify');

assert(isa(pgon_,'polygon'));

% test 'all'
pgon_ = compact(pgon,'all',0.1);

assert(isa(pgon_,'polygon'));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
