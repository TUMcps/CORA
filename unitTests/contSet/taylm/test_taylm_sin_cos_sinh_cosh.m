function res = test_taylm_sin_cos_sinh_cosh(~)
% test_taylm_sin_cos_sinh_cosh - unit-tests for Taylor models
%
% Syntax:  
%    res = test_taylm_sin_cos_sinh_cosh(~)
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
% Written:      15-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

    %% Test 1
    a = interval(0,2);
    a = taylm(a, 3); %-> 1 + x + [0,0]
    t = sin(a); %-> sin(1)*(1 - x^2/2) + cos(1)*(x - x^3/3!) + [0, 0.66667]
    eps = 10^-3;
    s = sin(1);
    c = cos(1);
    
    if ~appeq( getCoef(t),[s; c; -s/2; -c/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.04167), eps)
        res = 0;
        error('test 1 is failed')
    end   
    
    %% Test 2
    a = interval(0,2);
    a = taylm(a, 3); %-> 1 + x + [0,0]
    t = cos(a); %-> cos(1)*(1 - x^2/2) + sin(1)*(-x + x^3/3!) + [0, 0.66667]
    eps = 10^-3;
    s = sin(1);
    c = cos(1);
    
    if ~appeq( getCoef(t),[c; -s; -c/2; s/6], eps ) ||...
            ~appeq( getRem(t), interval(-0.01734, 0.04167), eps)
        res = 0;
        error('test 2 is failed')
    end
    
    %% Test 3
    a = interval(0,2);
    a = taylm(a, 3); %-> 1 + x + [0,0]
    t = sinh(a); %-> sinh(1) + cosh(1)*x + sinh(1)/2 *x^2 + cosh(1)/6 * x^4 + [0, 0.15112]
    eps = 10^-3;
    s = sinh(1);
    c = cosh(1);
    
    if ~appeq( getCoef(t),[s; c; s/2; c/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.15112), eps)
        res = 0;
        error('test 3 is failed')
    end 
    
    %% Test 4
    a = interval(0,2);
    a = taylm(a, 3); %-> 1 + x + [0,0]
    t = cosh(a); %-> cosh(1) + sinh(1)*x + cosh(1)/2 *x^2 + sinh(1)/6 * x^4 + [0, 0.15676]
    eps = 10^-3;
    s = sinh(1);
    c = cosh(1);
    
    if ~appeq( getCoef(t),[c; s; c/2; s/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.15676), eps)
        res = 0;
        error('test 4 is failed')
    end
    
    %% Test 5
    syms x
    a = taylm(1 + x,interval(-1,1), 3); %-> 1 + x + [0,0]
    t = sin(a); %-> sin(1)*(1 - x^2/2) + cos(1)*(x - x^3/3!) + [0, 0.66667]
    eps = 10^-3;
    s = sin(1);
    c = cos(1);
    
    if ~appeq( getCoef(t),[s; c; -s/2; -c/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.04167), eps)
        res = 0;
        error('test 5 is failed')
    end   
    
    %% Test 6
    syms x
    a = taylm(1 + x,interval(-1,1), 3); %-> 1 + x + [0,0]
    t = cos(a); %-> cos(1)*(1 - x^2/2) + sin(1)*(-x + x^3/3!) + [0, 0.66667]
    eps = 10^-3;
    s = sin(1);
    c = cos(1);
    
    if ~appeq( getCoef(t),[c; -s; -c/2; s/6], eps ) ||...
            ~appeq( getRem(t), interval(-0.01734, 0.04167), eps)
        res = 0;
        error('test 6 is failed')
    end
    
    %% Test 7
    syms x
    a = taylm(1 + x, interval(-1,1), 3); %-> 1 + x + [0,0]
    t = sinh(a); %-> sinh(1) + cosh(1)*x + sinh(1)/2 *x^2 + cosh(1)/6 * x^4 + [0, 0.15112]
    eps = 10^-3;
    s = sinh(1);
    c = cosh(1);
    
    if ~appeq( getCoef(t),[s; c; s/2; c/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.15112), eps)
        res = 0;
        error('test 7 is failed')
    end 
    
    %% Test 8
    syms x
    a = taylm(1 + x, interval(-1,1), 3); %-> 1 + x + [0,0]
    t = cosh(a); %-> cosh(1) + sinh(1)*x + cosh(1)/2 *x^2 + sinh(1)/6 * x^4 + [0, 0.15676]
    eps = 10^-3;
    s = sinh(1);
    c = cosh(1);
    
    if ~appeq( getCoef(t),[c; s; c/2; s/6], eps ) ||...
            ~appeq( getRem(t), interval(0, 0.15676), eps)
        res = 0;
        error('test 8 is failed')
    end  
    
    %% Test 9
    a = interval(0,pi/2);
    
    % loop over different maximum orders (sine)
    for i = 1:10
        t = taylm(a,i);
        int = interval(sin(t));
        if supremum(int) < 1 || infimum(int) > 0
           res = 0;
           error('test 9 is failed') 
        end
    end
    
    %% Test 10
    a = interval(0,pi/2);
    
    % loop over different maximum orders (cosine)
    for i = 1:10
        t = taylm(a,i);
        int = interval(cos(t));
        if supremum(int) < 1 || infimum(int) > 0
           res = 0;
           error('test 10 is failed') 
        end
    end
    
    
end

%------------- END OF CODE --------------