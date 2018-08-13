function res = test_taylm_asin_acos_atan(~)
% test_taylm_asin_acos_atan - unit-tests for Taylor models
%
% Syntax:  
%    res = test_taylm_asin_acos_atan(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none
%
% Author:       Dmitry Grebenyuk
% Written:      16-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

    %% Test 1
    a = interval(-1,1);
    a = taylm(a, 7); %-> x + [0,0]
    t = asin(a); %->  0 + x + x^3/3! + 3^2*x^5/5! + 3^2*5^2*x^7/7! + [0, 0.00092]
    eps = 10^-3;

    if ~appeq( getCoef(t),[0; 1; 1/6; 0.075; 0.0446], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.00092), eps)
        res = 0;
        error('test 1 is failed')
    end
    
    %% Test 2
    a = interval(-1,1);
    a = taylm(a, 7); %-> x + [0,0]
    t = acos(a); %->  pi/2 - x - x^3/3! - 3^2*x^5/5! - 3^2*5^2*x^7/7! + [-0.00092,0]
    eps = 10^-3;

    if ~appeq( getCoef(t), [pi/2; -1; -1/6; -0.075; -0.0446], eps ) ||...
            ~appeq( getRem(t), interval(-0.00092, 0), eps)
        res = 0;
        error('test 2 is failed')
    end
    
    %% Test 3
    a = interval(-1,1);
    a = taylm(a, 7); %-> x + [0,0]
    t = atan(a); %->  0 + x - x^3/3 + x^5/5 - x^7/7 + [-0.125,0.125]
    eps = 10^-3;

    if ~appeq( getCoef(t), [0; 1; -1/3; 1/5; -1/7], eps ) ||...
            ~appeq( getRem(t), interval(-0.125,0.125), eps)
        res = 0;
        error('test 3 is failed')
    end
    
    %% Test 4
    syms x
    a = taylm(x,interval(-1,1),7); %-> x + [0,0]
    t = asin(a); %->  0 + x + x^3/3! + 3^2*x^5/5! + 3^2*5^2*x^7/7! + [0, 0.00092]
    eps = 10^-3;

    if ~appeq( getCoef(t),[0; 1; 1/6; 0.075; 0.0446], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.00092), eps)
        res = 0;
        error('test 4 is failed')
    end
    
    %% Test 5
    syms x
    a = taylm(x,interval(-1,1), 7); %-> x + [0,0]
    t = acos(a); %->  pi/2 - x - x^3/3! - 3^2*x^5/5! - 3^2*5^2*x^7/7! + [-0.00092,0]
    eps = 10^-3;

    if ~appeq( getCoef(t), [pi/2; -1; -1/6; -0.075; -0.0446], eps ) ||...
            ~appeq( getRem(t), interval(-0.00092, 0), eps)
        res = 0;
        error('test 5 is failed')
    end
    
    %% Test 6
    syms x
    a = taylm(x, interval(-1,1), 7); %-> x + [0,0]
    t = atan(a); %->  0 + x - x^3/3 + x^5/5 - x^7/7 + [-0.125,0.125]
    eps = 10^-3;

    if ~appeq( getCoef(t), [0; 1; -1/3; 1/5; -1/7], eps ) ||...
            ~appeq( getRem(t), interval(-0.125,0.125), eps)
        res = 0;
        error('test 6 is failed')
    end
end

%------------- END OF CODE --------------