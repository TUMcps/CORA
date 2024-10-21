function res = test_finiteSignal_at
% test_finiteSignal_at - unit test of signal at function
%
% Syntax:
%    res = test_finiteSignal_at
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       09-December-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
sig1 = finiteSignal([1.2 2.3 3.0], [true false true]);
sig2 = finiteSignal(2.7, false);

% test
res = true;
assert(isequal(true, sig1.at(0)));
assert(isequal(true, sig1.at(1)));
assert(isequal(false, sig1.at(1.2)));
assert(isequal(false, sig1.at(2.2)));
assert(isequal(true, sig1.at(2.3)));
assert(isequal(true, sig1.at(3.0)));

assert(isequal(false, sig2.at(0)));
assert(isequal(false, sig2.at(1.9)));
assert(isequal(false, sig2.at(2.7)));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
