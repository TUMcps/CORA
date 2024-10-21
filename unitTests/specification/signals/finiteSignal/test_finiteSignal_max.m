function res = test_finiteSignal_max
% test_finiteSignal_max - unit test of signal max function
%
% Syntax:
%    res = test_finiteSignal_max
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
% Written:       23-January-2023
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
time1 = [1.2 2.4 2.7 3.9 4.0 4.5 6.7];
value1 = [tt uu ff uu ff tt uu];
sig1 = finiteSignal(time1, value1);

sig2 = finiteSignal(3.2, uu);

% test
assert(isequal(tt, max(sig1, interval(1.0, 2.0))));
assert(isequal(uu, max(sig1, interval(2.0, 2.2))));
assert(isequal(uu, max(sig1, interval(2.0, 3.0))));
assert(isequal(uu, max(sig1, interval(2.4, 2.6))));
assert(isequal(uu, max(sig1, interval(2.5, 2.7))));
assert(isequal(ff, max(sig1, interval(2.5, 2.6))));
assert(isequal(tt, max(sig1, interval(1.0, 5.0))));
assert(isequal(uu, max(sig1, interval(5.0, 9.0))));
assert(isequal(tt, max(sig1, interval(1.2, 2.0))));
assert(isequal(tt, max(sig1, interval(0.0, 6.7))));
assert(isequal(tt, max(sig1, interval(0.0, 9.0))));

assert(isequal(uu, max(sig2, interval(1.2, 2.3))));
assert(isequal(uu, max(sig2, interval(1.0, 1.0))));
assert(isequal(uu, max(sig2, interval(0.0, 3.2))));
assert(isequal(uu, max(sig2, interval(0.0, 5.0))));
assert(isequal(ff, max(sig2, interval(5.0, 9.0))));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
