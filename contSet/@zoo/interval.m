function int = interval(obj)
% interval - calculate the bounding interval of a zoo-object
%
% Syntax:  
%    int = interval(obj)
%
% Inputs:
%    obj - zoo object
%
% Outputs:
%    int - interval overapproximating the zoo (class interval)
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: zoo

% Author:       Niklas Kochdumper
% Written:      10-April-2018
% Last update:  --- 
% Last revision: ---

%------------- BEGIN CODE -------------

    int = arrayfun(@(a) s_zoo2int(a), obj, 'UniformOutput', 0);
    A = [int{:}];
    int = reshape(A, size(int)); 
    
end


%------------ Scalar function ---------------

function res = s_zoo2int(obj)

    res = interval(-inf,inf);
    for i = 1:length(obj.method)
       res = res & interval(obj.objects{i}); 
    end
    
end
    
%------------ END OF CODE ------------ 