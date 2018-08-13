function res = test_taylm_log(~)
% test_taylm_log - unit-test for Taylor models consisting of
%                       log functions
%
% Syntax:  
%    res = test_taylm_log(~)
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
% Written:      14-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

%% One-dimensional case
    a = interval(1,3);
    a = taylm(a, 3); %-> 2 + x + [0,0]
    t = log(a); %-> log(2) + x/2 - x^2/8 + x^3/24 + [-0.25, 0]
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [log(2); 1/2; -1/8; 1/24], eps ) ||...
            ~appeq( getRem(t), interval(-0.25,0), eps)
        res = 0;
        error('test 1 is failed')
    end 
    
    %% test 2
    syms x
    a = taylm(2 + x,interval(-1,1), 3); %-> 2 + x + [0,0]
    t = log(a); %-> log(2) + x/2 - x^2/8 + x^3/24 + [-0.25, 0]
    eps = 10^-3;
    
    if ~appeq( getCoef(t), [log(2); 1/2; -1/8; 1/24], eps ) ||...
            ~appeq( getRem(t), interval(-0.25,0), eps)
        res = 0;
        error('test 2 is failed')
    end  
end

%------------- END OF CODE --------------