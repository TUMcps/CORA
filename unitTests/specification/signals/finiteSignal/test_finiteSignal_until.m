function res = test_finiteSignal_until
% test_finiteSignal_until - unit test of signal until function
%
% Syntax:
%    res = test_finiteSignal_until
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
% Written:       17-August-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% signals
p = finiteSignal([1 3 4 14 16], [false true false true false]);
q = finiteSignal([4 6 12 16], [false true false true]);

% test
u1 = until(p, interval(2,4), q);
assert(isequal([8 12 16], u1.time));
assert(isequal([false true false], u1.value));

u2 = until(p, interval(1,2), q);
assert(isequal([4 5 10 13 16], u2.time));
assert(isequal([false true false true false], u2.value));


% signals
p = finiteSignal([2 11 14], [false true false]);
q = finiteSignal([8 10 14], [false true false]);

% test
u1 = until(p, interval(0,1), q);
assert(isequal([7 10 14], u1.time));
assert(isequal([false true false], u1.value));

u2 = until(p, interval(0,0), q);
assert(isequal([8 10 14], u2.time));
assert(isequal([false true false], u2.value));

u3 = until(p, interval(1,2), q);
assert(isequal([6 9 14], u3.time));
assert(isequal([false true false], u3.value));

u4 = until(p, interval(8,9), q);
assert(isequal(14, u4.time));
assert(isequal(false, u4.value));

u5 = until(p, interval(1,9), q);
assert(isequal([2 9 14], u5.time));
assert(isequal([false true false], u5.value));

% more tests
u6 = until(p, interval(6,7), q);
assert(isequal([2 4 14], u6.time));
assert(isequal([false true false], u6.value));

u7 = until(p, interval(0,4), q);
assert(isequal([4 10 14], u7.time));
assert(isequal([false true false], u7.value));

u8 = until(p, interval(1,4), q);
assert(isequal([4 9 14], u8.time));
assert(isequal([false true false], u8.value));

u9 = until(p, interval(0,5), q);
assert(isequal([3 10 14], u9.time));
assert(isequal([false true false], u9.value));


% signals
p = finiteSignal([2 11 14], [false true false]);
q = finiteSignal([1 4 8 10 14], [false true false true false]);

% test
u1 = until(p, interval(0,1), q);
assert(isequal([2 4 7 10 14], u1.time));
assert(isequal([false true false true false], u1.value));

u2 = until(p, interval(0,0), q);
assert(isequal([2 4 8 10 14], u2.time));
assert(isequal([false true false true false], u2.value));

u3 = until(p, interval(1,2), q);
assert(isequal([2 3 6 9 14], u3.time));
assert(isequal([false true false true false], u3.value));

u4 = until(p, interval(8,9), q);
assert(isequal(14, u4.time));
assert(isequal(false, u4.value));

u5 = until(p, interval(1,9), q);
assert(isequal([2 9 14], u5.time));
assert(isequal([false true false], u5.value));

% more tests
u6 = until(p, interval(6,7), q);
assert(isequal([2 4 14], u6.time));
assert(isequal([false true false], u6.value));

u7 = until(p, interval(0,4), q);
assert(isequal([2 10 14], u7.time));
assert(isequal([false true false], u7.value));

u8 = until(p, interval(1,4), q);
assert(isequal([2 3 4 9 14], u8.time));
assert(isequal([false true false true false], u8.value));

u9 = until(p, interval(0,5), q);
assert(isequal([2 10 14], u9.time));
assert(isequal([false true false], u9.value));


% signals
p = finiteSignal([1.5 8.5 10], [false true false]);
q = finiteSignal([2.0 6.5 7.5 9.0 10], [false true false true false]);

% test
u1 = until(p, interval(0,0.5), q);
assert(isequal([1.5 6.5 7.0 8.5 10], u1.time));
assert(isequal([false true false true false], u1.value));


% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
p = finiteSignal([3 5 14 17 19], [ff uu tt uu ff]);
q = finiteSignal([4 7 9 11 13 15 18 19], [ff uu tt uu tt ff uu ff]);

% test
u1 = until(p, interval(0,1), q, ff, uu);
assert(isequal([3 13 14 17 19], u1.time));
assert(isequal([ff uu ff uu ff], u1.value));

u2 = until(p, interval(0,1), q, ff, tt);
assert(isequal([6 9 10 13 19], u2.time));
assert(isequal([ff tt ff tt ff], u2.value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
