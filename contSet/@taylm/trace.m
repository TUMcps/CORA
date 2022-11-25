function res = trace(obj)
% trace - trace for TM matrices
%
% Syntax:  
%    res = trace(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      19-November-2017
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------

    res = obj(1,1);
    for i = 2:size(obj, 1)
        res = res + obj(i, i);
    end
    
end
%------------ END OF CODE ------------