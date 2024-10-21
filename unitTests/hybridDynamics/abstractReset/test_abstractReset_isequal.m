function res = test_abstractReset_isequal
% test_abstractReset_isequal - test function for equality check of
%    abstractReset objects
%
% Syntax:
%    res = test_abstractReset_isequal
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
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% standard case
reset1 = abstractReset(1,2,1);
reset2 = abstractReset(1,2,2);
assert(isequal(reset1,reset1));
assert(~isequal(reset1,reset2));
assert(~isequal(reset2,reset1));

% different operators
assert(reset1 == reset1);
assert(~(reset1 == reset2));
assert(~(reset2 == reset1));
assert(~(reset1 ~= reset1));
assert(reset1 ~= reset2);
assert(reset2 ~= reset1);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
