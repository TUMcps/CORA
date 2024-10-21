function res = test_emptySet_polygon
% test_emptySet_polygon - unit test function of polygon
%
% Syntax:
%    res = test_emptySet_polygon
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
% See also: -

% Authors:       Tobias Ladner
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
O = emptySet.empty(2);
pgon = polygon(O);
assert(representsa(pgon,"emptySet"));

% normal initialization
O = emptySet(2);
pgon = polygon(O);
assert(representsa(pgon,'emptySet'));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
