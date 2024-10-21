function res = test_polygon_mldivide()
% test_polygon_mldivide - unit test function for polygon/mldivide
%
% Syntax:
%    res = test_polygon_mldivide()
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
% Written:       11-October-2024
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

% check hole
pgon_I = polygon(interval([0.4;0.4],[0.75;0.75]));
pgon_hole = pgon \ pgon_I;

% compute vertices
V_true = [ ...
 0.9527821496566151, 0.9102352155402841, 0.8929320288328670, 0.8407431981130701, 0.5226947790160240, 0.4428188422351334, 0.1995708604681839, 0.0595160380235712, 0.2651894190237225, 0.3123440066941408, 0.0222097785726014, 0.0344695245434839, 0.3030947130117874, 0.4086261330703079, 0.7687159962242812, 0.9538774735922314, 0.8985957224020654,              NaN, 0.4000000000000000, 0.7500000000000000, 0.7500000000000000, 0.4000000000000000 ; ...
 0.4492815327156741, 0.3280519505009173, 0.0452718951436371, 0.0159086201489517, 0.0243094359597330, 0.0015051689556644, 0.0736017167258955, 0.4230814965883369, 0.5222028944090500, 0.7706254812339077, 0.9104892159525940, 0.9509129441148924, 0.9284666248925266, 0.9842482902762220, 0.9476667564545247, 0.9640790006827448, 0.5779175416460717,              NaN, 0.4000000000000000, 0.4000000000000000, 0.7500000000000000, 0.7500000000000000 ; ...
 ];
V = vertices(pgon_hole);

% filter nan values
V(isnan(V)) = [];
V_true(isnan(V_true)) = [];
assert(compareMatrices(V_true,V));

% check with empty polytope
pgon_hole = pgon \ polygon.empty(2);
assert(isequal(pgon,pgon_hole));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
