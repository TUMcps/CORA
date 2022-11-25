function res = atan( obj )
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
% See also: taylm, acos, asin
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      16-August-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_atan(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_atan( obj )
    
    if obj.monomials(1,1) == 0                 
        c_f = obj.coefficients(1);
    else
        c_f = obj.coefficients( obj.monomials(:,1) == 0 );
    end
   
    % producing T_tilda from the manual
    if isempty(c_f)
        T = obj;
        c_f = 0;
    else    
        T = (obj - c_f) / (1 + c_f .* obj);
    end
    
    % factorial initialisation
    factor1 = -1;
    T_factor = 1;
    % sum initialisation
    T1 = atan(c_f);

    for i = 1:obj.max_order
        T_factor = T_factor .* T;
        if mod(i,2) == 1
            factor1 = factor1 .* -1;
            T1 = T1 + T_factor .* factor1 ./ i;
        end
    end   

    rem = interval(T);
    
    % temporaly increase max-order to get tight bounds for interval "remPow"
    T_factor.max_order = T_factor.max_order + T.max_order;
    remPow = interval(T_factor.*T);     % remPow = B( T^(k+1) )
    
    res = T1;
    
    powr = (obj.max_order + 1);
    
    res.remainder = res.remainder +...
        1 ./ powr * remPow .* ...
        cos(atan(interval(0,1) .* rem)).^powr .*...
        sin(powr .* (atan(interval(0,1) .* rem) + pi/2));
end
%------------- END OF CODE --------------

