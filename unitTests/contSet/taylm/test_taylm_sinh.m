function res = test_taylm_sinh
% test_taylm_sinh - unit test of hyperbolic sine function
%
% Syntax:
%    res = test_taylm_sinh
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
t = sinh(a); %-> sinh(1) + cosh(1)*x + sinh(1)/2 *x^2 + cosh(1)/6 * x^4 + [0, 0.15112]
eps = 10^-3;
s = sinh(1);
c = cosh(1);

if ~appeq( getCoef(t),[s; c; s/2; c/6], eps ) ||...
        ~appeq( getRem(t), interval(0, 0.15112), eps)
    res = false;
end 

%% Test 2
syms x
a = taylm(1 + x, interval(-1,1), 3); %-> 1 + x + [0,0]
t = sinh(a); %-> sinh(1) + cosh(1)*x + sinh(1)/2 *x^2 + cosh(1)/6 * x^4 + [0, 0.15112]
eps = 10^-3;
s = sinh(1);
c = cosh(1);

if ~appeq( getCoef(t),[s; c; s/2; c/6], eps ) ||...
        ~appeq( getRem(t), interval(0, 0.15112), eps)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
