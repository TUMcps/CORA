function res = test_taylm_exp(~)
% test_taylm_exp - unit-tests for Taylor models consisting of
%                       exponent functions
%
% Syntax:  
%    res = test_taylm_exp(~)
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

%% One-dimensional case
    a = interval(0,2);
    a = taylm(a, 3); %-> 1 + x + [0,0]
    t = exp(a); %-> exp(1)*(1 + x + x^2/2 + x^3/3!) + [0, 0.30788]
    eps = 10^-3;
    
    if ~appeq( getCoef(t), exp(1)*[1; 1; 1/2; 1/6], eps ) ||...
            ~appeq( getRem(t), interval(0,0.30788), eps)
        res = 0;
        error('test 1 is failed')
    end   
    
    %% test 2
    syms x
    a = taylm(1 + x,interval(-1,1), 3); %-> 1 + x + [0,0]
    t = exp(a); %-> exp(1)*(1 + x + x^2/2 + x^3/3!) + [0, 0.30788]
    eps = 10^-3;
    
    if ~appeq( getCoef(t), exp(1)*[1; 1; 1/2; 1/6], eps ) ||...
            ~appeq( getRem(t), interval(0,0.30788), eps)
        res = 0;
        error('test 2 is failed')
    end 
end

%------------- END OF CODE --------------