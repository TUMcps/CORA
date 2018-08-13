function res = tan( obj )
% tan - Overloaded 'tan()' operator for a Taylor model
%
% Syntax:  
%    res = tan(obj)
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
% See also: taylm, sin, cos
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      15-August-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_tan(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_tan( obj )

    if isempty(obj.monomials)    % Taylor model without polynomial part

        res = obj;
        res.remainder = tan(obj.remainder);

    else                    % Taylor model with polynomial part
        res = sin(obj) ./ cos(obj);
    end
end

%------------- END OF CODE --------------

