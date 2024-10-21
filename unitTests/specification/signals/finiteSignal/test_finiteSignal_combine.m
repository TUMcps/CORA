function res = test_finiteSignal_combine
% test_finiteSignal_combine - unit test of signal combine function
%
% Syntax:
%    res = test_finiteSignal_combine
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

time2 = [0.9 1.3 3.0 3.9 4.2 6.3];
value2 = [uu ff tt uu tt uu];
sig2 = finiteSignal(time2, value2);

sigF = finiteSignal(5, ff);
sigT = finiteSignal(5, tt);
sigU = finiteSignal(5, uu);

% conjunction
timeCon = [0.9 1.3 2.4 2.7 3.9 4.0 4.2 6.3];
valueCon = [uu ff uu ff uu ff tt uu];

con = sig1 & sig2;

assert(isequal(timeCon, con.time));
assert(isequal(valueCon, con.value));

assert(isequal(sigF, sigF & sig1));
assert(isequal(sigF, sigF & sig2));

% disjunction
timeDis = [1.2 1.3 3.0 3.9 4.5 6.3];
valueDis = [tt uu tt uu tt uu];

dis = sig1 | sig2;

assert(isequal(timeDis, dis.time));
assert(isequal(valueDis, dis.value));

assert(isequal(sigT, sigT | sig1));
assert(isequal(sigT, sigT | sig2));

% negation
neg1 = ~ sig1;
neg2 = ~ sig2;

assert(isequal(sig1.time, neg1.time));
assert(isequal([ff uu tt uu tt ff uu], neg1.value));

assert(isequal(sig2.time, neg2.time));
assert(isequal([uu tt ff uu ff uu], neg2.value));

assert(isequal(sigF, ~ sigT));
assert(isequal(sigT, ~ sigF));
assert(isequal(sigU, ~ sigU));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
