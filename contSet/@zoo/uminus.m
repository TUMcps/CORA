function obj = uminus(obj)
% uminus - Overloaded '-' operator for single operand
%
% Syntax:  
%    res = uplus(obj)
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

% Author:       Dmitry Grebenyuk
% Written:      06-November-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    obj = arrayfun(@(a) s_uminus(a), obj, 'UniformOutput', 0);
    A = cat(1, obj{:});
    obj = reshape(A, size(obj));

end

%% --------------- Implementation for a scalar --------------

function obj = s_uminus(obj)

    for i = 1:length(obj.method)
       obj.objects{i} = -obj.objects{i}; 
    end
end
%------------- END OF CODE --------------