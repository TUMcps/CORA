function res = test_taylm_sqrt
% test_taylm_sqrt - unit test of square root function
%
% Syntax:
%    res = test_taylm_sqrt
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
% Written:       14-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

%% One-dimensional case
a = interval(1,3);
a = taylm(a, 4); %-> 2 + x + [0,0]
t = sqrt(a); %-> sqrt(2)(1 + x/4 - x^2/32 + 3*x^3/384 - 15*x^4/6144) + [-0.02734, 0.02734]
eps = 10^-3;

if ~appeq( getCoef(t), sqrt(2)*[1; 1/4; -1/32; 3/384; -15/6144], eps ) ||...
        ~appeq( getRem(t), interval(-0.02734, 0.02734), eps)
    res = false;
end  

%% test 2
syms x
a = taylm(2 + x,interval(-1,1), 4); %-> 2 + x + [0,0]
t = sqrt(a); %-> sqrt(2)(1 + x/4 - x^2/32 + 3*x^3/384 - 15*x^4/6144) + [-0.02734, 0.02734]
eps = 10^-3;

if ~appeq( getCoef(t), sqrt(2)*[1; 1/4; -1/32; 3/384; -15/6144], eps ) ||...
        ~appeq( getRem(t), interval(-0.02734, 0.02734), eps)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
