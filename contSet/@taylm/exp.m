function res = exp(obj)
% exp - Overloaded 'exp()' operator for a Taylor model
%
% Syntax:  
%    res = exp(obj)
%
% Inputs:
%    obj - (input value) a taylm object
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      31-July-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_exp(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_exp(obj)
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
        T = obj - c_f;
    end
    % sum initialisation
    T1 = 1;
    % factorial initialisation
    factor = 1;
    T_factor = 1;
    
    for i = 1:obj.max_order
        factor = factor * i;
        T_factor = T_factor .* T;
        T1 = T1 + T_factor ./ factor;
    end
    
    excf = exp(c_f);
    rem = interval(T);
    
    % temporaly increase max-order to get tight bounds for interval "remPow"
    T_factor.max_order = T_factor.max_order + T.max_order;
    remPow = interval(T_factor.*T);     % remPow = B( T^(k+1) )
    
    res = T1 * excf;
    res.remainder = res.remainder +...
        excf /(factor * (obj.max_order + 1)) *...
        ( remPow .* exp( interval(0,1) .* rem ) );
end

%------------- END OF CODE --------------