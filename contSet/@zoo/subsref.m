function res = subsref(obj, S)
% subsref - Overloads the opertor that selects elements, e.g. T(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    res = subsref(obj, S)
%
% Inputs:
%    obj - a zoo object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    res - element or elemets of the taylm matrix
%
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: zoo

% Author:       Niklas Kochdumper
% Written:      10-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = obj;
    
    % loop over all combined calls to subsref
    for j = 1:length(S)
        
        % override subsref for type obj(i,j)
        if strcmp(S(j).type,'()')
            
            % loop over all objects stored in the zoo
            for i = 1:length(res.method)
                res.objects{i} = subsref(res.objects{i}, S(j));
            end
            
        % for all other types, call the build in function
        else
            res = builtin('subsref', res, S(j));
        end
    end
end
%------------- END OF CODE --------------