function res = test_taylm_cos
% test_taylm_cos - unit test of cosine function
%
% Syntax:
%    res = test_taylm_cos
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
t = cos(a); %-> cos(1)*(1 - x^2/2) + sin(1)*(-x + x^3/3!) + [0, 0.66667]
eps = 10^-3;
s = sin(1);
c = cos(1);

if ~appeq( getCoef(t),[c; -s; -c/2; s/6], eps ) ||...
        ~appeq( getRem(t), interval(-0.01734, 0.04167), eps)
    res = false;
end

%% Test 2
syms x
a = taylm(1 + x,interval(-1,1), 3); %-> 1 + x + [0,0]
t = cos(a); %-> cos(1)*(1 - x^2/2) + sin(1)*(-x + x^3/3!) + [0, 0.66667]
eps = 10^-3;
s = sin(1);
c = cos(1);

if ~appeq( getCoef(t),[c; -s; -c/2; s/6], eps ) ||...
        ~appeq( getRem(t), interval(-0.01734, 0.04167), eps)
    res = false;
end

%% Test 3
a = interval(0,pi/2);

% loop over different maximum orders (cosine)
for i = 1:10
    t = taylm(a,i);
    int = interval(cos(t));
    if supremum(int) < 1 || infimum(int) > 0
       res = false;
    end
end

% ------------------------------ END OF CODE ------------------------------
