function res = test_matZonotope_printSet
% test_matZonotope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_matZonotope_printSet
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

% test normal set
C = [0 0; 0 0];
G(:,:,1) = [1 3; -1 2];
G(:,:,2) = [2 0; 1 -1];
matZ = matZonotope(C,G);

printSet(matZ);
printSet(matZ,'high');
printSet(matZ,'high',true);
printSet(matZ,'high',false);

% test fid
filename = 'test.txt';
printSet(filename,matZ,'high',true);
matZ_copy = eval(fileread(filename));
assert(isequal(matZ,matZ_copy));
delete(filename)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
