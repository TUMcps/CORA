function res = test_signal_cutoff
% test_signal_cutoff - unit test of signal cutoff function
%
% Syntax:
%    res = test_signal_cutoff
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
% Written:       19-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% signals
sig1 = signal([1.2 2.3 3.0], [true false true]);
sig2 = signal(2.7, false);

% test
c1 = cutoff(sig1, 2.7);

res(end+1,1) = isequal([1.2 2.3 2.7], c1.time);
res(end+1,1) = isequal([true false true], c1.value);

c2 = cutoff(sig1, 3.0);

res(end+1,1) = isequal([1.2 2.3 3.0], c2.time);
res(end+1,1) = isequal([true false true], c2.value);

c3 = cutoff(sig1, 2.0);

res(end+1,1) = isequal([1.2 2.0], c3.time);
res(end+1,1) = isequal([true false], c3.value);

c4 = cutoff(sig1, 2.3);

res(end+1,1) = isequal([1.2 2.3], c4.time);
res(end+1,1) = isequal([true false], c4.value);

c5 = cutoff(sig1, 0.1);

res(end+1,1) = isequal(0.1, c5.time);
res(end+1,1) = isequal(true, c5.value);

c6 = cutoff(sig1, 0);

res(end+1,1) = isequal(0, c6.time);
res(end+1,1) = isequal(true, c6.value);

c7 = cutoff(sig2, 2.7);

res(end+1,1) = isequal(2.7, c7.time);
res(end+1,1) = isequal(false, c7.value);

c8 = cutoff(sig2, 1.2);

res(end+1,1) = isequal(1.2, c8.time);
res(end+1,1) = isequal(false, c8.value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
