function obj = uminus(obj)
% uminus - Overloaded '-' operator for single operand
%
% Syntax:
%    res = uminus(obj)
%    res = -obj
%
% Inputs:
%    obj - a taylm object
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

% Authors:       Dmitry Grebenyuk
% Written:       15-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    obj = arrayfun(@(a) aux_s_uminus(a), obj, 'UniformOutput', 0);
    A = cat(1, obj{:});
    obj = reshape(A, size(obj));

end


% Auxiliary functions -----------------------------------------------------

% Implementation for a scalar
function obj = aux_s_uminus(obj)

    obj.coefficients = -obj.coefficients;
    obj.remainder = -obj.remainder;

end

% ------------------------------ END OF CODE ------------------------------
