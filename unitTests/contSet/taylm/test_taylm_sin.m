function res = test_taylm_sin
% test_taylm_sin - unit test of sine function
%
% Syntax:
%    res = test_taylm_sin
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none

% Authors:       Dmitry Grebenyuk
% Written:       15-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

%% Test 1
a = interval(0,2);
a = taylm(a, 3); %-> 1 + x + [0,0]
t = sin(a); %-> sin(1)*(1 - x^2/2) + cos(1)*(x - x^3/3!) + [0, 0.66667]
eps = 10^-3;
s = sin(1);
c = cos(1);

assert(appeq( getCoef(t),[s; c; -s/2; -c/6], eps ))
assert(appeq( getRem(t), interval(0, 0.04167), eps))


%% Test 2
syms x
a = taylm(1 + x,interval(-1,1), 3); %-> 1 + x + [0,0]
t = sin(a); %-> sin(1)*(1 - x^2/2) + cos(1)*(x - x^3/3!) + [0, 0.66667]
eps = 10^-3;
s = sin(1);
c = cos(1);

assert(appeq( getCoef(t),[s; c; -s/2; -c/6], eps ))
assert(appeq( getRem(t), interval(0, 0.04167), eps))


%% Test 3
a = interval(0,pi/2);

% loop over different maximum orders (sine)
for i = 1:10
    t = taylm(a,i);
    int = interval(sin(t));
    assertLoop(supremum(int) >= 1,i)
    assertLoop(infimum(int) <= 0,i)
end

% ------------------------------ END OF CODE ------------------------------
