function res = power(base,exponent)
% power - Overloaded '.^' operator for taylm (power)
%
% Syntax:  
%    res = power(base,exponent)
%
% Inputs:
%    base - taylm object
%    exponent - taylm object
%
% Outputs:
%    res - taylm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      06-August-2017
% Last update:  ---               
% Last revision:---

%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_power(a,exponent), base, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_power(base,exponent)

    res = base;

    % trivial case
    if exponent == 0
        res = 1;
        
    elseif exponent > 0   
        
        % positive integer exponent
        if mod(exponent,1) == 0
            for i = 1:(exponent-1)
                res = res .* base;
            end
            
        % multiple of 1/2 => sqrt()
        elseif mod(2*exponent,1) == 0
            base = sqrt(base);
            res = base;
            for i = 1:(2*exponent-1)
                res = res .* base;
            end
            
        % positive real exponent
        else
            error('Operation not implemented for class "taylm"!');
        end
        
    % negative exponent
    else
        res = inverse(base).^ -exponent;
    end

end

%------------ END OF CODE ------------