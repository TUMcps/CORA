function res = testLong_taylm_division
% testLong_taylm_division - unit-tests for Taylor models consisting
%    of division operation
%
% Syntax:
%    res = testLong_taylm_division
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
% Written:       06-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% loop over different intervals
for i = 1:5
    
    % create a random interval 
    temp = rand(2,1)*i;
    D = interval(min(temp),max(temp));
    
    for j = 1:20
        % create the taylor model
        t = taylm(D,j,'x');
        T = repmat(t,1,j);
        
        % create a taylor polynomial
        e = 1:j;
        coeff = rand(j,1);
        temp = ((T.^e) * coeff);
        poly = temp(1);
        for k = 2:length(temp)
           poly = poly + temp(k); 
        end
        
        % make sure that the taylor polynomial does not contain 0
        temp = interval(poly);
        if contains(temp,0)
           r = 2*rand();
           if r > 1
              poly = poly - 1.1*infimum(temp);  
           else
              poly = poly - 1.1*supremum(temp); 
           end
        end
       
        % determine the bounds using Taylor models
        int = interval(1/poly);
        
        % determine the bounds using random sampling
        x = -1:0.01:1;
        y = zeros(size(x));
        
        for k = 1:length(poly.coefficients)
            y = y + poly.coefficients(k)*x.^poly.monomials(k,1);
        end
        intReal = 1/interval(min(y),max(y));
        
        % compare the results
        if infimum(int) > infimum(intReal) || supremum(int) < supremum(intReal)
            throw(CORAerror('CORA:testFailed'));
        end          
    end
end

% ------------------------------ END OF CODE ------------------------------
