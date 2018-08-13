function [newObj] = subsref(obj, S)
% subsref - Overloads the opertor that selects elements, e.g. I(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    [element] = subsref(obj, S)
%
% Inputs:
%    obj - interval object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    element - element or elemets of the interval matrix
%
% Example: 
%    a=interval([-1 1], [1 2]);
%    a(1,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-June-2015 
% Last update:  22-June-2015 
% Last revision:---

%------------- BEGIN CODE --------------

%obtain sub-intervals from the interval object
newObj = obj;
%check if parantheses are used to select elements
if strcmp(S.type,'()')
    % only one index specified
    if length(S.subs)==1
        newObj.inf=obj.inf(S.subs{1});
        newObj.sup=obj.sup(S.subs{1});
    %two indices specified
    elseif length(S.subs)==2
        %Select column of obj
        if strcmp(S.subs{1},':')
            column=S.subs{2};
            newObj.inf=obj.inf(:,column);
            newObj.sup=obj.sup(:,column);
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row=S.subs{1};
            newObj.inf=obj.inf(row,:);
            newObj.sup=obj.sup(row,:);
        %Select single element of V    
        elseif isnumeric(S.subs{1}) & isnumeric(S.subs{1})
            row=S.subs{1};
            column=S.subs{2};
            newObj.inf=obj.inf(row,column);
            newObj.sup=obj.sup(row,column);
        end
    end
end

%------------- END OF CODE --------------