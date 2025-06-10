function res = test_zonotope_printSet
% test_zonotope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_zonotope_printSet
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
Z = zonotope.empty(2);

printSet(Z);
printSet(Z,'high');
printSet(Z,'high',true);
printSet(Z,'high',false);

% test normal set
c = [1;1];
G = [1 1 1; 1 -1 0];
Z = zonotope(c,G);

printSet(Z);
printSet(Z,'high');
printSet(Z,'high',true);
printSet(Z,'high',false);

% test fid
filename = 'test.txt';
printSet(filename,Z,'high',true);
Z_copy = eval(fileread(filename));
assert(isequal(Z,Z_copy));
delete(filename)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
