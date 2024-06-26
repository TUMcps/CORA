function obj = uminus(obj)
% uminus - Overloaded '-' operator for single operand
%
% Syntax:
%    res = uminus(obj)
%
% Inputs:
%    obj - a zoo object
%
% Outputs:
%    res - a zoo object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, interval

% Authors:       Dmitry Grebenyuk
% Written:       06-November-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    obj = arrayfun(@(a) aux_s_uminus(a), obj, 'UniformOutput', false);
    A = cat(1, obj{:});
    obj = reshape(A, size(obj));

end


% Auxiliary functions -----------------------------------------------------

% Implementation for a scalar

function obj = aux_s_uminus(obj)

    for i = 1:length(obj.method)
       obj.objects{i} = -obj.objects{i}; 
    end
end

% ------------------------------ END OF CODE ------------------------------
