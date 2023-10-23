function val = subsref(intMat,S)
% subsref - Overloads the operator that selects elements, e.g. A(1,2),
%    where the element of the first row and second column is referred to.
%
% Syntax:
%    element = subsref(intMat,S)
%
% Inputs:
%    intMat - intervalMatrix object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    val - element or elemets of the interval hull matrix
%
% Example:
%    C = [0 2; 3 1];
%    D = [1 2; 1 1];
%    intMat = intervalMatrix(C,D);
%    intMat(2,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       23-September-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%check if parantheses are used to select elements
if strcmp(S.type,'()')
    if length(S.subs)==2
        %Select column of V
        if strcmp(S.subs{1},':')
            column=S.subs{2};
            val=intMat.int(:,column);
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row=S.subs{1};
            val=intMat.int(row,:);
        %Select single element of V    
        elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
            row=S.subs{1};
            column=S.subs{2};
            val=intMat.int(row,column);
        end
    %no selection if elements not proper specified  
    else
        val=[];
    end
    
%check if dot is used to select elements
elseif strcmp(S.type,'.')
    if isprop(intMat,S.subs)
        val = intMat.(S.subs);
    else
        throw(CORAerror('CORA:specialError',...
            strcat('Object of class "',class(intMat),'" does not have a property "',S.subs,'".')));
    end
%no selection if parantheses are not used    
else
    val=[];
end

% ------------------------------ END OF CODE ------------------------------
