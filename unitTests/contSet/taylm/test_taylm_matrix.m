function res = test_taylm_matrix(~)
% test_taylm_matrix - unit-tests for Taylor models
%
% Syntax:  
%    res = test_taylm_matrix(~)
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
% Written:      15-October-2017
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

    %% Test 1
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = eye(2)*[a b; b a];
    eps = 10^-6;
    
    if ~appeq( getCoef(t(1,1)),[2; 1], eps ) ||...
            ~appeq( getCoef(t(1,2)),[-2; 1], eps ) ||...
            ~appeq( getCoef(t(2,1)),[-2; 1], eps ) ||...
            ~appeq( getCoef(t(2,2)),[2; 1], eps ) ||...
            ~appeq( getRem(t), interval(0, 0), eps)
        res = 0;
        error('test 1 is failed')
    end
    
    %% Test 2
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = eye(2).*[a b; b a];
    eps = 10^-6;
    
    if ~appeq( getCoef(t(1,1)),[2; 1], eps ) ||...
            ~appeq( getCoef(t(1,2)),[0], eps ) ||...
            ~appeq( getCoef(t(2,1)),[0], eps ) ||...
            ~appeq( getCoef(t(2,2)),[2; 1], eps ) ||...
            ~appeq( getRem(t), interval(0, 0), eps)
        res = 0;
        error('test 2 is failed')
    end
    
    %% Test 3
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = [a b; b a]^0;
    
    if t ~= eye(2)
        res = 0;
        error('test 3 is failed')
    end
    
    %% Test 4
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = [a b; b a]^0;
    eps = 10^-6;
    
    if t ~= ones(2,2)
        res = 0;
        error('test 4 is failed')
    end
    
     %% Test 5
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = [a b; b a]*[a b; b a];
    t11 = a*a + b*b;
    t12 = a*b + b*a;
    t21 = b*a + a*b;
    t22 = b*b + a*a;
    eps = 10^-6;
    
    if ~appeq( getCoef(t(1,1)),getCoef(t11), eps ) ||...
            ~appeq( getCoef(t(1,2)),getCoef(t12), eps ) ||...
            ~appeq( getCoef(t(2,1)),getCoef(t21), eps ) ||...
            ~appeq( getCoef(t(2,2)),getCoef(t22), eps ) ||...
            ~appeq( getRem(t), interval(0, 0), eps)
        res = 0;
        error('test 5 is failed')
    end
    
     %% Test 6
    a = taylm(interval(1, 3), 7, {'a'}); %-> a + 2 + [0,0]
    b = taylm(interval(-3, -1), 7, {'b'}); %-> b - 2 + [0,0]
    t = [a b; b a].*[a b; b a];
    t11 = a.*a ;
    t12 = b.*b;
    t21 = b.*b;
    t22 = a.*a;
    eps = 10^-6;
    
    if ~appeq( getCoef(t(1,1)),getCoef(t11), eps ) ||...
            ~appeq( getCoef(t(1,2)),getCoef(t12), eps ) ||...
            ~appeq( getCoef(t(2,1)),getCoef(t21), eps ) ||...
            ~appeq( getCoef(t(2,2)),getCoef(t22), eps ) ||...
            ~appeq( getRem(t), interval(0, 0), eps)
        res = 0;
        error('test 6 is failed')
    end
    
     %% Test 7
    syms x1 x2
    A = [1 + 0.1*x1, 1; 1, -2 - 0.5*x2];
    B = taylm(A,interval([-1;-1],[1;1]), 6);
    C = (1/B);
    t = C*B;
    
    eps = 10^-6;
    
    if ~appeq( getCoef(t(1,1)),1, eps ) ||...
            ~appeq( getCoef(t(1,2)),0, eps ) ||...
            ~appeq( getCoef(t(2,1)),0, eps ) ||...
            ~appeq( getCoef(t(2,2)),1, eps )
        res = 0;
        error('test 7 is failed')
    end
    
end

%------------- END OF CODE --------------