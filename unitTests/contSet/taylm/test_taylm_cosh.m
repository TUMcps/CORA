function res = test_taylm_cosh
% test_taylm_cosh - unit tests of hyperbolic cosine function
%
% Syntax:
%    res = test_taylm_cosh
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
t = cosh(a); %-> cosh(1) + sinh(1)*x + cosh(1)/2 *x^2 + sinh(1)/6 * x^4 + [0, 0.15676]
eps = 10^-3;
s = sinh(1);
c = cosh(1);

if ~appeq( getCoef(t),[c; s; c/2; s/6], eps ) ||...
        ~appeq( getRem(t), interval(0, 0.15676), eps)
    res = false;
end

%% Test 2
syms x
a = taylm(1 + x, interval(-1,1), 3); %-> 1 + x + [0,0]
t = cosh(a); %-> cosh(1) + sinh(1)*x + cosh(1)/2 *x^2 + sinh(1)/6 * x^4 + [0, 0.15676]
eps = 10^-3;
s = sinh(1);
c = cosh(1);

if ~appeq( getCoef(t),[c; s; c/2; s/6], eps ) ||...
        ~appeq( getRem(t), interval(0, 0.15676), eps)
    res = false;
end  

% ------------------------------ END OF CODE ------------------------------
