function res = asin( obj )
% asin - Overloaded 'asin()' operator for a taylm expression
%
% Syntax:  
%    res = asin(obj)
%
% Inputs:
%    obj - a taylm expression
%
% Outputs:
%    res - a taylm expression object
%
% Example: 
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, acos
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      16-August-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_asin(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_asin( obj )
    if obj.monomials(1,1) == 0                 
        c_f = obj.coefficients(1);
    else
        c_f = obj.coefficients( obj.monomials(:,1) == 0 );
    end
   
    % check if the input taylor model is admissible (> 0)
    temp = interval(obj);
    if infimum(temp) < -1 || supremum(temp) > 1
       error('Function ''taylm/asin'' is not defined for values <-1 and >1!'); 
    end
    
    % producing T_tilda from the manual
    if isempty(c_f)
        T = obj;
        c_f = 0;
    else    
        T = obj .* sqrt(1 - c_f.^2) - c_f .* sqrt(1 - obj.^2);
    end
    
    % factorial initialisation
    factor = 1;
    factor1 = 1;
    T_factor = 1;
    % sum initialisation
    T1 = asin(c_f);

    for i = 1:obj.max_order
        factor = factor * i;
        T_factor = T_factor .* T;
        if mod(i,2) == 1
            factor1 = factor1 .* dasin(i);
            T1 = T1 + T_factor .* factor1 ./ factor;
        end
    end   

    rem = interval(T);
    
    % temporaly increase max-order to get tight bounds for interval "remPow"
    T_factor.max_order = T_factor.max_order + T.max_order;
    remPow = interval(T_factor.*T);     % remPow = B( T^(k+1) )
    
    res = T1;
    
    factor = factor * (obj.max_order + 1);
    
    res.remainder = res.remainder +...
        1 /factor * remPow .*...
        asin(interval(0,1) .* rem).^(obj.max_order + 1);
end

function res = dasin(k)
    if k >= 3
        res = (k - 2)^2;
    else
        res = 1;
    end
end
%------------- END OF CODE --------------
