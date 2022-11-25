function res = test_taylm_atan
% test_taylm_atan - unit test of inverse tangent function
%
% Syntax:  
%    res = test_taylm_atan
%
% Inputs:
%    -
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
    t = atan(a); %->  0 + x - x^3/3 + x^5/5 - x^7/7 + [-0.125,0.125]
    eps = 10^-3;

    if ~appeq( getCoef(t), [0; 1; -1/3; 1/5; -1/7], eps ) ||...
            ~appeq( getRem(t), interval(-0.125,0.125), eps)
        res = 0;
        error('test 1 is failed')
    end
    
    %% Test 2
    syms x
    a = taylm(x, interval(-1,1), 7); %-> x + [0,0]
    t = atan(a); %->  0 + x - x^3/3 + x^5/5 - x^7/7 + [-0.125,0.125]
    eps = 10^-3;

    if ~appeq( getCoef(t), [0; 1; -1/3; 1/5; -1/7], eps ) ||...
            ~appeq( getRem(t), interval(-0.125,0.125), eps)
        res = 0;
        error('test 2 is failed')
    end
end

%------------- END OF CODE --------------