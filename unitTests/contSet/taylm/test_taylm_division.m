function res = test_taylm_division(~)
% test_taylm_division - unit-tests for Taylor models consisting of 
%                             division operation
%
% Syntax:  
%    res = test_taylm_divison(~)
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
% Written:      06-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

%% One-dimensional case
    %% Test 1
    a = interval(1,3);
    t = taylm(a, 3); %-> 2 + x + [0,0]
    inv_t = 1/t; %-> 1/2 - x/4 + x^2/8 - x^3/16 + [0,1]
    eps = 10^-3;

    if ~appeq(getCoef(inv_t), [0.5; -0.25; 0.125; -0.0625], eps) ||...
            ~appeq(getRem(inv_t), interval(0,1), eps)
        res = 0;
        error('test 1 is failed')
    end

    %% Test 2
    a = interval(1,3);
    t = taylm(a, 3) + a; %-> 2 + x + [1,3]
    t = t/3; %-> 2/3 + x/3 + [1/3,1]
    eps = 10^-3;

    if ~appeq(getCoef(t), [2/3; 1/3], eps) ||...
            ~appeq(getRem(t), interval(1/3,1), eps)
        res = 0;
        error('test 2 is failed')
    end
    
    %% Test 3
    syms x
    t = taylm(2 + x,interval(-1,1), 3); %-> 2 + x + [0,0]
    inv_t = 1/t; %-> 1/2 - x/4 + x^2/8 - x^3/16 + [0,1]
    eps = 10^-3;

    if ~appeq(getCoef(inv_t), [0.5; -0.25; 0.125; -0.0625], eps) ||...
            ~appeq(getRem(inv_t), interval(0,1), eps)
        res = 0;
        error('test 3 is failed')
    end

    %% Test 4
    syms x
    t = taylm(2 + x, interval(1, 3), 3); %-> 4 + x + [0,0]
    t = t/3; %-> 4/3 + x/3 + [0,0]
    eps = 10^-3;

    if ~appeq(getCoef(t), [4/3; 1/3], eps) ||...
            ~appeq(getRem(t), interval(0,0), eps)
        res = 0;
        error('test 4 is failed')
    end
    
    %% Test 5: Random tests
    
    tol = 1e-7;
    
    % loop over different intervals
    for i = 1:5
        
        % create a random interval 
        temp = rand(2,1)*i;
        iTemp = interval(min(temp),max(temp));
        
        for j = 1:20
            % create the taylor model
            t = taylm(iTemp,j);
            T = repmat(t,1,j);
            
            % create a taylor polynomial
            e = 1:j;
            coeff = rand(j,1);
            temp = ((T.^e) * coeff);
            poly = temp(1);
            for k = 2:length(temp)
               poly = poly + temp(k); 
            end
           
            % determine the bounds
            int = interval(1/poly);
            intReal = 1/interval(poly);
            
            if infimum(intReal) - infimum(int) < -tol || supremum(int) - supremum(intReal) < -tol
                res = 0;
                error('test 5 is failed');
            end          
        end
    end
end

%------------- END OF CODE --------------