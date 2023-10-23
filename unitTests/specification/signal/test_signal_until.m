function res = test_signal_until
% test_signal_until - unit test of signal until function
%
% Syntax:
%    res = test_signal_until
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% signals
p = signal([1 3 4 14 16], [false true false true false]);
q = signal([4 6 12 16], [false true false true]);

% test
u1 = until(p, interval(2,4), q);
res(end+1,1) = isequal([8 12 16], u1.time);
res(end+1,1) = isequal([false true false], u1.value);

u2 = until(p, interval(1,2), q);
res(end+1,1) = isequal([4 5 10 13 16], u2.time);
res(end+1,1) = isequal([false true false true false], u2.value);


% signals
p = signal([2 11 14], [false true false]);
q = signal([8 10 14], [false true false]);

% test
u1 = until(p, interval(0,1), q);
res(end+1,1) = isequal([7 10 14], u1.time);
res(end+1,1) = isequal([false true false], u1.value);

u2 = until(p, interval(0,0), q);
res(end+1,1) = isequal([8 10 14], u2.time);
res(end+1,1) = isequal([false true false], u2.value);

u3 = until(p, interval(1,2), q);
res(end+1,1) = isequal([6 9 14], u3.time);
res(end+1,1) = isequal([false true false], u3.value);

u4 = until(p, interval(8,9), q);
res(end+1,1) = isequal(14, u4.time);
res(end+1,1) = isequal(false, u4.value);

u5 = until(p, interval(1,9), q);
res(end+1,1) = isequal([2 9 14], u5.time);
res(end+1,1) = isequal([false true false], u5.value);

u6 = until(p, interval(6,7), q);
res(end+1,1) = isequal([2 4 14], u6.time);
res(end+1,1) = isequal([false true false], u6.value);

u7 = until(p, interval(0,4), q);
res(end+1,1) = isequal([4 10 14], u7.time);
res(end+1,1) = isequal([false true false], u7.value);

u8 = until(p, interval(1,4), q);
res(end+1,1) = isequal([4 9 14], u8.time);
res(end+1,1) = isequal([false true false], u8.value);

u9 = until(p, interval(0,5), q);
res(end+1,1) = isequal([3 10 14], u9.time);
res(end+1,1) = isequal([false true false], u9.value);


% signals
p = signal([2 11 14], [false true false]);
q = signal([1 4 8 10 14], [false true false true false]);

% test
u1 = until(p, interval(0,1), q);
res(end+1,1) = isequal([2 4 7 10 14], u1.time);
res(end+1,1) = isequal([false true false true false], u1.value);

u2 = until(p, interval(0,0), q);
res(end+1,1) = isequal([2 4 8 10 14], u2.time);
res(end+1,1) = isequal([false true false true false], u2.value);

u3 = until(p, interval(1,2), q);
res(end+1,1) = isequal([2 3 6 9 14], u3.time);
res(end+1,1) = isequal([false true false true false], u3.value);

u4 = until(p, interval(8,9), q);
res(end+1,1) = isequal(14, u4.time);
res(end+1,1) = isequal(false, u4.value);

u5 = until(p, interval(1,9), q);
res(end+1,1) = isequal([2 9 14], u5.time);
res(end+1,1) = isequal([false true false], u5.value);

u6 = until(p, interval(6,7), q);
res(end+1,1) = isequal([2 4 14], u6.time);
res(end+1,1) = isequal([false true false], u6.value);

u7 = until(p, interval(0,4), q);
res(end+1,1) = isequal([2 10 14], u7.time);
res(end+1,1) = isequal([false true false], u7.value);

u8 = until(p, interval(1,4), q);
res(end+1,1) = isequal([2 3 4 9 14], u8.time);
res(end+1,1) = isequal([false true false true false], u8.value);

u9 = until(p, interval(0,5), q);
res(end+1,1) = isequal([2 10 14], u9.time);
res(end+1,1) = isequal([false true false], u9.value);


% signals
p = signal([1.5 8.5 10], [false true false]);
q = signal([2.0 6.5 7.5 9.0 10], [false true false true false]);

% test
u1 = until(p, interval(0,0.5), q);
res(end+1,1) = isequal([1.5 6.5 7.0 8.5 10], u1.time);
res(end+1,1) = isequal([false true false true false], u1.value);


% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
p = signal([3 5 14 17 19], [ff uu tt uu ff]);
q = signal([4 7 9 11 13 15 18 19], [ff uu tt uu tt ff uu ff]);

% test
u1 = until(p, interval(0,1), q, ff, uu);
res(end+1,1) = isequal([3 13 14 17 19], u1.time);
res(end+1,1) = isequal([ff uu ff uu ff], u1.value);

u2 = until(p, interval(0,1), q, ff, tt);
res(end+1,1) = isequal([6 9 10 13 19], u2.time);
res(end+1,1) = isequal([ff tt ff tt ff], u2.value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
