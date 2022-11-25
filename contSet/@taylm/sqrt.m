function res = sqrt( obj )
% sqrt - compute formula the square root for Taylor models
%
% Syntax:  
%    res = sqrt(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      14-August-2017
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_sqrt(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_sqrt( obj )        
    if isempty(obj.monomials)    % Taylor model without polynomial part
        
        res = obj;
        res.remainder = sqrt(obj.remainder);
        
    else                    % Taylor model with polynomial part
        
        if obj.monomials(1, 1) == 0                 
            c_f = obj.coefficients(1);
        else
            c_f = obj.coefficients( obj.monomials(:,1) == 0 );
        end
        
        %producing T_tilda from the manual
        if isempty(c_f)
            T = obj;
            c_f = 0;
        else    
            T = obj - c_f;
        end
        
        % check if the input taylor model is admissible (> 0)
        rem = interval(T);
        temp = rem + c_f;
        if infimum(temp) < 0
           error('Function ''taylm/sqrt'' is not defined for negative inputs!'); 
        end
    
        % sum initialisation
        T1 = 1;
        sqrt_c_f = sqrt(c_f);
        c_f = 1./c_f;
        factor = 1;
        factor1 = 1;
        factor2 = 1;
        T_factor = 1;

        for i = 1:obj.max_order
            factor = factor .* i;
            factor1 = factor1 .* 2;
            factor2 = factor2 .* two_k_m_3(i);
            T_factor = T_factor .* T .* c_f; 
            T1 = T1 + (-1).^(i-1) .* T_factor .* factor2 ./ factor1 ./ factor;
        end
        
        factor = factor .* (obj.max_order + 1);
        factor1 = factor1 .* 2;
        factor2 = factor2 .* two_k_m_3(obj.max_order + 1);
        
        % temporaly increase max-order to get tight bounds for interval "remPow"
        T_factor.max_order = T_factor.max_order + T.max_order;
        remPow = interval(T_factor.*T.*c_f);     % remPow = B( (T/c_f)^(k+1) )

        res = sqrt_c_f .* T1;
    	res.remainder = res.remainder +...
            (-1)^(obj.max_order) .* sqrt_c_f .* factor2 ./ factor1 ./ factor .*...
            remPow ./ (1 + interval(0,1) .* rem .* c_f).^(obj.max_order + 1/2);
        
    end
end

function res = two_k_m_3(k)
    if k > 1
        res = (2 .* k - 3);
    else
        res = 1;
    end
end

%------------ END OF CODE ------------


