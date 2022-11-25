function res = tanh( obj )
% tanh - Overloaded 'tanh()' operator for a Taylor model
%
% Syntax:  
%    res = tanh(obj)
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
% See also: taylm, sinh, cosh
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      04-September-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_tanh(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_tanh( obj )
    if isempty(obj.monomials)    % Taylor model without polynomial part

        res = obj;
        res.remainder = tanh(obj.remainder);

    else                    % Taylor model with polynomial part
        res = sinh(obj) ./ cosh(obj);
    end
end

%------------- END OF CODE --------------
