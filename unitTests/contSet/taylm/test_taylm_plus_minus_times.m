function res = test_taylm_plus_minus_times(~)
% test_taylm_plus_minus_times - unit-tests for Taylor models consisting of 
%                             plus, munus, and times operations
%
% Syntax:  
%    res = test_taylm_plus_minus_times(~)
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
% Written:      07-August-2017
%               14-Octover-2017 (DG) test for syms initialization are added
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

    %% Test 1
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3); %-> 0.5 + 0.5x + [0,0]
    y = taylm(y, 3); %-> 0.5 + 0.5y + [0,0]
    t = x - x*y; %-> 1/4 + 1/4*x - 1/4*y - 1/4*x*y  + [0, 0]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; -1/4; 1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(0,0), eps)
        res = 0;
        error('test 1 is failed')
    end
    
    %% Test 2
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3); %-> 0.5 + 0.5x + [0,0]
    y = taylm(y, 3); %-> 0.5 + 0.5y + [0,0]
    t = x + x*y; %-> 2/4 + 3/4*x + 1/4*y + 1/4*x*y  + [0, 0]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [3/4; 1/4; 3/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(0,0), eps)
        res = 0;
        error('test 2 is failed')
    end
    
    %% Test 3
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = x + x*y; %-> 3/4 + 3/4*x + 1/4*y + 1/4*x*y  + [-13, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [3/4; 1/4; 3/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12,4), eps)
        res = 0;
        error('test 3 is failed')
    end
    
    %% Test 4
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = x - x*y; %-> 3/4 + 3/4*x + 1/4*y + 1/4*x*y  + [-13, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; -1/4; 1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(-2,14), eps)
        res = 0;
        error('test 4 is failed')
    end
    
    %% Test 5
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = 2*x*y; %-> 2/4 + 2/4*x + 2/4*y + 2/4*x*y  + [-24, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [2/4; 2/4; 2/4; 2/4], eps ) ||...
            ~appeq( getRem(t), interval(-24,4), eps)
        res = 0;
        error('test 5 is failed')
    end
    
    %% Test 6
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = x*y - 2; %-> -7/4 + 2/4*x + 2/4*y + 2/4*x*y  + [-12, 2]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [-7/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12,2), eps)
        res = 0;
        error('test 6 is failed')
    end
    
    %% Test 7
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = 2 - x*y; %-> 7/4 - 1/4*x - 1/4*y - 1/4*x*y  + [-2, 12]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [7/4; -1/4; -1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(-2,12), eps)
        res = 0;
        error('test 7 is failed')
    end
    
    %% Test 8
    x = interval(0,1);
    y = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    t = 2*y*x - x*y; %-> 1/4 + 1/4*x + 1/4*y + 1/4*x*y  + [-26, 16]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-26, 16), eps)
        res = 0;
        error('test 8 is failed')
    end
    
    %% Test 9
    x = interval(0,1);
    y = interval(0,1);
    z = interval(0,1);
    x = taylm(x, 3) + 2*x; %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y; %-> 1/2 + 1/2*y + [-4,0]
    z = taylm(z, 3); %-> 1/2 + 1/2*z + [0,0]
    t = 2*y*x*z - x*y; %-> 0 + 1/4*z + 1/4*y*z + 1/4*x*z + 1/4*x*y*z + [-26, 16]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [0; 1/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-26, 16), eps)
        res = 0;
        error('test 9 is failed')
    end
    
    %% Test 10
    x = interval(0,1);
    y = interval(0,1);
    a = interval(0,1);
    b = interval(0,1);
    x = taylm(x, 3) + 2*x;    %-> 1/2 + 1/2*x + [0,2]
    y = taylm(y, 3) - 4*y;    %-> 1/2 + 1/2*y + [-4,0]
    a = taylm(a, 3);          %-> 1/2 + 1/2*a + [0,0]
    b = taylm(b, 3);          %-> 1/2 + 1/2*b + [0,0]
    t = x*y - a*b;   %-> 0 + 1/4*x + 1/4*y + 1/4*x*y 
                    %     + 1/4*a + 1/4*b + 1/4*a*b + [-12, 2]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [0; -1/4; -1/4; 1/4; 1/4; -1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12, 2), eps)
        res = 0;
        error('test 10 is failed')
    end
    
    %% Test 11
    x = interval(-2,0);
    x = -taylm(x, 2); %-> 1 - x + [0,0]
    t = x*x*x; %-> 1 - 3x + 3x^2  + [-1, 1]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1; -3; 3], eps ) ||...
            ~appeq( getRem(t), interval(-1,1), eps)
        res = 0;
        error('test 11 is failed')
    end
    
     %% Test 12
    syms x y
    p1 = 0.5 + 0.5*x;
    p2 = 0.5 + 0.5*y;
    x = taylm(p1,interval(-1,1), 3); %-> 0.5 + 0.5x + [0,0]
    y = taylm(p2,interval(-1,1), 3); %-> 0.5 + 0.5y + [0,0]
    t = x - x*y; %-> 1/4 + 1/4*x - 1/4*y - 1/4*x*y  + [0, 0]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; -1/4; 1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(0,0), eps)
        res = 0;
        error('test 12 is failed')
    end
    
    %% Test 13
    syms x y
    p1 = 0.5 + 0.5*x;
    p2 = 0.5 + 0.5*y;
    x = taylm(p1,interval(-1,1), 3); %-> 0.5 + 0.5x + [0,0]
    y = taylm(p2,interval(-1,1), 3); %-> 0.5 + 0.5y + [0,0]
    t = x + x*y; %-> 2/4 + 3/4*x + 1/4*y + 1/4*x*y  + [0, 0]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [3/4; 1/4; 3/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(0,0), eps)
        res = 0;
        error('test 13 is failed')
    end
    
    %% Test 14
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = x + x*y; %-> 3/4 + 3/4*x + 1/4*y + 1/4*x*y  + [-13, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [3/4; 1/4; 3/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12,4), eps)
        res = 0;
        error('test 14 is failed')
    end
    
    %% Test 15
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = x - x*y; %-> 3/4 + 3/4*x + 1/4*y + 1/4*x*y  + [-13, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; -1/4; 1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(-2,14), eps)
        res = 0;
        error('test 15 is failed')
    end
    
    %% Test 16
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = 2*x*y; %-> 2/4 + 2/4*x + 2/4*y + 2/4*x*y  + [-24, 4]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [2/4; 2/4; 2/4; 2/4], eps ) ||...
            ~appeq( getRem(t), interval(-24,4), eps)
        res = 0;
        error('test 16 is failed')
    end
    
    %% Test 17
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = x*y - 2; %-> -7/4 + 2/4*x + 2/4*y + 2/4*x*y  + [-12, 2]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [-7/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12,2), eps)
        res = 0;
        error('test 17 is failed')
    end
    
    %% Test 18
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = 2 - x*y; %-> 7/4 - 1/4*x - 1/4*y - 1/4*x*y  + [-2, 12]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [7/4; -1/4; -1/4; -1/4], eps ) ||...
            ~appeq( getRem(t), interval(-2,12), eps)
        res = 0;
        error('test 18 is failed')
    end
    
    %% Test 19
    syms x y
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    t = 2*y*x - x*y; %-> 1/4 + 1/4*x + 1/4*y + 1/4*x*y  + [-26, 16]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-26, 16), eps)
        res = 0;
        error('test 19 is failed')
    end
    
    %% Test 20
    syms x y z
    p1 = 1/2 + 1/2*x;
    p2 = 1/2 + 1/2*y;
    p3 = 1/2 + 1/2*z;
    x = taylm(p1, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(p2, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    z = taylm(p3, interval(-1,1), 3); %-> 1/2 + 1/2*z + [0,0]
    t = 2*y*x*z - x*y; %-> 0 + 1/4*z + 1/4*y*z + 1/4*x*z + 1/4*x*y*z + [-26, 16]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [0; 1/4; 1/4; 1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-26, 16), eps)
        res = 0;
        error('test 20 is failed')
    end
    
    %% Test 21
    syms x y a b
    x = taylm(1/2 + 1/2*x, interval(-1,1), 3) + interval(0,2); %-> 1/2 + 1/2*x + [0,2]
    y = taylm(1/2 + 1/2*y, interval(-1,1), 3) + interval(-4,0); %-> 1/2 + 1/2*y + [-4,0]
    a = taylm(1/2 + 1/2*a, interval(-1,1), 3);          %-> 1/2 + 1/2*a + [0,0]
    b = taylm(1/2 + 1/2*b, interval(-1,1), 3);          %-> 1/2 + 1/2*b + [0,0]
    t = x*y - a*b;   %-> 0 + 1/4*x + 1/4*y + 1/4*x*y 
                    %     + 1/4*a + 1/4*b + 1/4*a*b + [-12, 2]

    eps = 10^-3;
    
    if ~appeq( getCoef(t), [0; -1/4; -1/4; 1/4; 1/4; -1/4; 1/4], eps ) ||...
            ~appeq( getRem(t), interval(-12, 2), eps)
        res = 0;
        error('test 21 is failed')
    end
    
    %% Test 22
    syms x
    x = -taylm(-1 + x,interval(-1,1), 2); %-> 1 - x + [0,0]
    t = x*x*x; %-> 1 - 3x + 3x^2  + [-1, 1]
    
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [1; -3; 3], eps ) ||...
            ~appeq( getRem(t), interval(-1,1), eps)
        res = 0;
        error('test 22 is failed')
    end
end

%------------- END OF CODE -----------