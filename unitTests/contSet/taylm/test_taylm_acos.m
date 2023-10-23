function res = test_taylm_acos
% test_taylm_acos - unit test of inverse cosine function
%
% Syntax:
%    res = test_taylm_acos
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
% Written:       16-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
    
%% Test 1
a = interval(-1,1);
a = taylm(a, 7); %-> x + [0,0]
t = acos(a); %->  pi/2 - x - x^3/3! - 3^2*x^5/5! - 3^2*5^2*x^7/7! + [-0.00092,0]
eps = 10^-3;

if ~appeq( getCoef(t), [pi/2; -1; -1/6; -0.075; -0.0446], eps ) ||...
        ~appeq( getRem(t), interval(-0.00092, 0), eps)
    res = false;
end

%% Test 2
syms x
a = taylm(x,interval(-1,1), 7); %-> x + [0,0]
t = acos(a); %->  pi/2 - x - x^3/3! - 3^2*x^5/5! - 3^2*5^2*x^7/7! + [-0.00092,0]
eps = 10^-3;

if ~appeq( getCoef(t), [pi/2; -1; -1/6; -0.075; -0.0446], eps ) ||...
        ~appeq( getRem(t), interval(-0.00092, 0), eps)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
