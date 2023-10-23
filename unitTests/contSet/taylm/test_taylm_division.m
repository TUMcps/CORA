function res = test_taylm_division
% test_taylm_division - unit-tests for Taylor models consisting of 
%    division operation
%
% Syntax:
%    res = test_taylm_division
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
% Written:       06-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
%% One-dimensional case
%% Test 1
a = interval(1,3);
t = taylm(a, 3); %-> 2 + x + [0,0]
inv_t = 1/t; %-> 1/2 - x/4 + x^2/8 - x^3/16 + [0,1]
eps = 10^-3;

if ~appeq(getCoef(inv_t), [0.5; -0.25; 0.125; -0.0625], eps) ||...
        ~appeq(getRem(inv_t), interval(0,1), eps)
    res = false;
end

%% Test 2
a = interval(1,3);
t = taylm(a, 3) + a; %-> 2 + x + [1,3]
t = t/3; %-> 2/3 + x/3 + [1/3,1]
eps = 10^-3;

if ~appeq(getCoef(t), [2/3; 1/3], eps) ||...
        ~appeq(getRem(t), interval(1/3,1), eps)
    res = false;
end

%% Test 3
syms x
t = taylm(2 + x,interval(-1,1), 3); %-> 2 + x + [0,0]
inv_t = 1/t; %-> 1/2 - x/4 + x^2/8 - x^3/16 + [0,1]
eps = 10^-3;

if ~appeq(getCoef(inv_t), [0.5; -0.25; 0.125; -0.0625], eps) ||...
        ~appeq(getRem(inv_t), interval(0,1), eps)
    res = false;
end

%% Test 4
syms x
t = taylm(2 + x, interval(1, 3), 3); %-> 4 + x + [0,0]
t = t/3; %-> 4/3 + x/3 + [0,0]
eps = 10^-3;

if ~appeq(getCoef(t), [4/3; 1/3], eps) ||...
        ~appeq(getRem(t), interval(0,0), eps)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
