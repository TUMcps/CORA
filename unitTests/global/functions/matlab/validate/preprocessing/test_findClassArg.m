function res = test_findClassArg
% test_findClassArg - unit test function for finding class arguments
%
% Syntax:
%    res = test_findClassArg()
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

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate class objects
Z = zonotope([2;-1],[2 4 1; 0 -2 1]);
I = interval([3;-1],[4;2]);

[obj1,obj2] = findClassArg(Z,I,'zonotope');
res = obj1 == Z && obj2 == I;
[obj1,obj2] = findClassArg(Z,I,'interval');
res(end+1,1) = obj1 == I && obj2 == Z;
[obj1,obj2] = findClassArg(I,Z,'interval');
res(end+1,1) = obj1 == I && obj2 == Z;
try
    [obj1,obj2] = findClassArg(Z,I,'polyZonotope');
    res(end+1,1) = false;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
