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

    int = arrayfun(@(a) aux_zoo2int(a), obj, 'UniformOutput', false);
    A = [int{:}];
    int = reshape(A, size(int)); 
    
end


%------------ Scalar function ---------------

function res = aux_zoo2int(obj)

    res = interval(-inf,inf);
    for i = 1:length(obj.method)
       res = and_(res,interval(obj.objects{i},'exact')); 
    end
    
end
    
%------------ END OF CODE ------------ 