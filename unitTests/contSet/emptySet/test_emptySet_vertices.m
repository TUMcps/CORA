function res = test_emptySet_vertices
% test_emptySet_vertices - unit test function of vertices
%
% Syntax:
%    res = test_emptySet_vertices
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
% Written:       22-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% init empty set
n = 2;
O = emptySet(n);

% check vertices size
V = vertices(O);
assert(isequal(size(V),[2,0]))

end

% ------------------------------ END OF CODE ------------------------------
