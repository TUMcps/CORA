function res = test_stlInterval_isequal
% test_stlInterval_isequal - unit test function of isequal
%
% Syntax:
%    res = test_stlInterval_isequal
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

% Authors:       Florian Lercher
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty interval
I1 = stlInterval();
I2 = stlInterval(0,1);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));

% closed/open interval
I1 = stlInterval(0,1);
I2 = stlInterval(0,1,false);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));

% unbounded interval
I1 = stlInterval(0,inf);
I2 = stlInterval(0,42);
assert(isequal(I1,I1));
assert(~isequal(I1,I2));

res = true;

% ------------------------------ END OF CODE ------------------------------
