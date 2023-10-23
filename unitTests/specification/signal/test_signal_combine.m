function res = test_signal_combine
% test_signal_combine - unit test of signal combine function
%
% Syntax:
%    res = test_signal_combine
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
% Written:       11-August-2022
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

sigF = signal(5, ff);
sigT = signal(5, tt);
sigU = signal(5, uu);

% conjunction
timeCon = [0.9 1.3 2.4 2.7 3.9 4.0 4.2 6.3];
valueCon = [uu ff uu ff uu ff tt uu];

con = sig1 & sig2;

res(end+1,1) = isequal(timeCon, con.time);
res(end+1,1) = isequal(valueCon, con.value);

res(end+1,1) = isequal(sigF, sigF & sig1);
res(end+1,1) = isequal(sigF, sigF & sig2);

% disjunction
timeDis = [1.2 1.3 3.0 3.9 4.5 6.3];
valueDis = [tt uu tt uu tt uu];

dis = sig1 | sig2;

res(end+1,1) = isequal(timeDis, dis.time);
res(end+1,1) = isequal(valueDis, dis.value);

res(end+1,1) = isequal(sigT, sigT | sig1);
res(end+1,1) = isequal(sigT, sigT | sig2);

% negation
neg1 = ~ sig1;
neg2 = ~ sig2;

res(end+1,1) = isequal(sig1.time, neg1.time);
res(end+1,1) = isequal([ff uu tt uu tt ff uu], neg1.value);

res(end+1,1) = isequal(sig2.time, neg2.time);
res(end+1,1) = isequal([uu tt ff uu ff uu], neg2.value);

res(end+1,1) = isequal(sigF, ~ sigT);
res(end+1,1) = isequal(sigT, ~ sigF);
res(end+1,1) = isequal(sigU, ~ sigU);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
