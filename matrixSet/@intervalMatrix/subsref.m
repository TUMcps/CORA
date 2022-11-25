function [element] = subsref(obj, S)
% subsref - Overloads the operator that selects elements, e.g. A(1,2),
%    where the element of the first row and second column is referred to.
%
% Syntax:  
%    [element] = subsref(obj, S)
%
% Inputs:
%    obj - interval matrix object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    element - element or elemets of the interval hull matrix
%
% Example: 
%    A=intervalMatrix(center,delta);
%    IA(2,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      23-September-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%check if parantheses are used to select elements
if strcmp(S.type,'()')
    if length(S.subs)==2
        %Select column of V
        if strcmp(S.subs{1},':')
            column=S.subs{2};
            element=obj.int(:,column);
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row=S.subs{1};
            element=obj.int(row,:);
        %Select single element of V    
        elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
            row=S.subs{1};
            column=S.subs{2};
            element=obj.int(row,column);
        end
    %no selection if elements not proper specified  
    else
        element=[];
    end
    
%check if dot is used to select elements
elseif strcmp(S.type,'.')
    if isprop(obj,S.subs)
        element = obj.(S.subs);
    else
        error(strcat('Object of class "',class(obj),'" does not have a property "',S.subs,'".'));
    end
%no selection if parantheses are not used    
else
    element=[];
end

%------------- END OF CODE --------------