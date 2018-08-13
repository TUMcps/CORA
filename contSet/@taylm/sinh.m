function res = sinh( obj )
% sinh - Overloaded 'sinh()' operator for a Taylor model
%
% Syntax:  
%    res = sinh(obj)
%
% Inputs:
%    obj - a taylm object
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, cosh
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      15-August-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_sinh(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_sinh( obj )
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
    
    % factorial initialisation
    factor = 1;
    T_factor = 1;
    cs_cf = cosh(c_f);
    sn_cf = sinh(c_f);
    % sum initialisation
    T1 = sn_cf;

    for i = 1:obj.max_order
        factor = factor * i;
        T_factor = T_factor .* T;
        if mod(i,2) == 0
            T1 = T1 + sn_cf .* T_factor ./ factor;
        else
            T1 = T1 + cs_cf .* T_factor ./ factor;
        end
    end   

    rem = interval(T);
    
    % temporaly increase max-order to get tight bounds for interval "remPow"
    T_factor.max_order = T_factor.max_order + T.max_order;
    remPow = interval(T_factor.*T);     % remPow = B( T^(k+1) )
    
    res = T1;
    
    if mod(obj.max_order + 1, 2) == 0
        J0 = sinh(c_f + interval(0,1) .* rem);
    else
        J0 = cosh(c_f + interval(0,1) .* rem);
    end
    
    res.remainder = res.remainder +...
        1 /(factor * (obj.max_order + 1)) *...
        remPow .* J0;
end

%------------- END OF CODE --------------