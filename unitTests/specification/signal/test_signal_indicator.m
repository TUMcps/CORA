function res = test_signal_indicator
% test_signal_indicator - unit test of signal indicator function
%
% Syntax:
%    res = test_signal_indicator
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
% Written:       16-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
sig1 = signal.indicator(3.0, interval(2.3, 2.7), tt);
sig2 = signal.indicator(3.0, interval(2.3, 2.7), uu);
sig3 = signal.indicator(3.0, interval(2.3, 2.7), ff);

sig4 = signal.indicator(3.0, interval(0.0, 2.7), tt);
sig5 = signal.indicator(3.0, interval(2.3, 3.0), tt);
sig6 = signal.indicator(3.0, interval(0.0, 3.0), tt);

sig7 = signal.indicator(3.0, interval(2.7, 4.3), tt);
sig8 = signal.indicator(3.0, interval(0, 4.3), tt);
sig9 = signal.indicator(3.0, interval(-1.3, 1.4), tt);

sig10 = signal.indicator(3.0, interval(4.2, 5.7), tt);
sig11 = signal.indicator(3.0, interval(-2.7, -2.4), tt);
sig12 = signal.indicator(3.0, interval(-2.3, 3.7), tt);

% test
res(end+1,1) = isequal([2.3 2.7 3.0], sig1.time);
res(end+1,1) = isequal([ff tt ff], sig1.value);

res(end+1,1) = isequal([2.3, 2.7, 3.0], sig2.time);
res(end+1,1) = isequal([ff uu ff], sig2.value);

res(end+1,1) = isequal(3.0, sig3.time);
res(end+1,1) = isequal(ff, sig3.value);

res(end+1,1) = isequal([2.7 3.0], sig4.time);
res(end+1,1) = isequal([tt ff], sig4.value);

res(end+1,1) = isequal([2.3 3.0], sig5.time);
res(end+1,1) = isequal([ff tt], sig5.value);

res(end+1,1) = isequal(3.0, sig6.time);
res(end+1,1) = isequal(tt, sig6.value);

res(end+1,1) = isequal([2.7 3.0], sig7.time);
res(end+1,1) = isequal([ff tt], sig7.value);

res(end+1,1) = isequal(3.0, sig8.time);
res(end+1,1) = isequal(tt, sig8.value);

res(end+1,1) = isequal([1.4 3.0], sig9.time);
res(end+1,1) = isequal([tt ff], sig9.value);

res(end+1,1) = isequal(3.0, sig10.time);
res(end+1,1) = isequal(ff, sig10.value);

res(end+1,1) = isequal(3.0, sig11.time);
res(end+1,1) = isequal(ff, sig11.value);

res(end+1,1) = isequal(3.0, sig12.time);
res(end+1,1) = isequal(tt, sig12.value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
