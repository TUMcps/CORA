function res = test_signal_findIntervals
% test_signal_findIntervals - unit test of signal findIntervals function
%
% Syntax:
%    res = test_signal_findIntervals
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

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
time1 = [1.2 2.4 2.7 3.9 4.0 4.5 6.7];
value1 = [tt uu ff uu ff tt uu];
sig1 = signal(time1, value1);

time2 = [0.9 1.3 3.0 3.9 4.2 6.3];
value2 = [uu ff tt uu tt uu];
sig2 = signal(time2, value2);

sig3 = signal(3.2, uu);

% test
int1 = findIntervals(sig1, @(v) v > ff);
res(end+1,1) = isequal(3, length(int1));
res(end+1,1) = isequal(interval(0, 2.4), int1(1));
res(end+1,1) = isequal(interval(2.7, 3.9), int1(2));
res(end+1,1) = isequal(interval(4.0, 6.7), int1(3));

int2 = findIntervals(sig1, @(v) v == tt);
res(end+1,1) = isequal(2, length(int2));
res(end+1,1) = isequal(interval(0, 1.2), int2(1));
res(end+1,1) = isequal(interval(4.0, 4.5), int2(2));

int3 = findIntervals(sig1, @(v) true);
res(end+1,1) = isequal(1, length(int3));
res(end+1,1) = isequal(interval(0, 6.7), int3(1));

int4 = findIntervals(sig1, @(v) false);
res(end+1,1) = isequal(0, length(int4));

int5 = findIntervals(sig2, @(v) v > ff);
res(end+1,1) = isequal(2, length(int5));
res(end+1,1) = isequal(interval(0, 0.9), int5(1));
res(end+1,1) = isequal(interval(1.3, 6.3), int5(2));

int6 = findIntervals(sig2, @(v) v == tt);
res(end+1,1) = isequal(2, length(int6));
res(end+1,1) = isequal(interval(1.3, 3.0), int6(1));
res(end+1,1) = isequal(interval(3.9, 4.2), int6(2));

int7 = findIntervals(sig3, @(v) v > ff);
res(end+1,1) = isequal(1, length(int7));
res(end+1,1) = isequal(interval(0, 3.2), int7(1));

int8 = findIntervals(sig3, @(v) v == tt);
res(end+1,1) = isequal(0, length(int8));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
