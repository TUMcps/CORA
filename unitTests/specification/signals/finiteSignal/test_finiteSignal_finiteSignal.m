function res = test_finiteSignal_finiteSignal
% test_finiteSignal_finiteSignal - unit test function for constructor
%
% Syntax:
%    res = test_finiteSignal_finiteSignal
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

% Authors:       Benedikt Seidl
% Written:       12-May-2023
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

assertThrowsAs(@finiteSignal,'MATLAB:validators:mustBeNonnegative',[-1 3 4],[true false true]);

assertThrowsAs(@finiteSignal,'CORA:wrongInputInConstructor',[],true);

assertThrowsAs(@finiteSignal,'CORA:wrongInputInConstructor',1,[]);

assertThrowsAs(@finiteSignal,'CORA:wrongInputInConstructor',[1 2 3],[true false]);

assertThrowsAs(@finiteSignal,'CORA:wrongInputInConstructor',[3 2 1],[true false true]);

assertThrowsAs(@finiteSignal,'CORA:wrongInputInConstructor',[1 2 2 3],[true false true false]);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
