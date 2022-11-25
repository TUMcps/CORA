function res = log( obj )
% log - compute formula 'log' for Taylor models
%
% Syntax:  
%    res = log(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval, interval, taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      13-August-2017
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_log(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_log( obj )
        
    if isempty(obj.monomials)    % Taylor model without polynomial part
        
        res = obj;
        res.remainder = log(obj.remainder);
        
    else                    % Taylor model with polynomial part
        
        if obj.monomials(1,1) == 0                 
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
           error('Function ''taylm/log'' is not defined for values <= 0!'); 
        end
    
        % sum initialisation
        T1 = log(c_f);
        c_f = 1./c_f;
        T_factor = 1;

        for i = 1:obj.max_order
            T_factor = T_factor .* T .* c_f; 
            T1 = T1 + (-1).^(i+1) .* T_factor ./ i;
        end

        
        % temporaly increase max-order to get tight bounds for interval "remPow"
        T_factor.max_order = T_factor.max_order + T.max_order;
        remPow = interval(T_factor.*T.*c_f);     % remPow = B( (T/c_f)^(k+1) )
        
        res = T1;
    	res.remainder = res.remainder +...
            (-1)^(obj.max_order + 2) ./ (obj.max_order + 1) .* remPow ./ ...
            (1 + interval(0,1) .* rem .* c_f).^(obj.max_order + 1);
        
    end
end

%------------ END OF CODE ------------
