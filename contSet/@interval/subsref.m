function newObj = subsref(I,S)
% subsref - Overloads the operator that selects elements, e.g. I(1,2),
%    where the element of the first row and second column is referred to.
%
% Syntax:
%    newObj = subsref(I,S)
%
% Inputs:
%    I - interval object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    newObj - element or elemets of the interval matrix
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

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       19-June-2015 
% Last update:   22-June-2015 
%                12-November-2018 (NK, default to build in for other cases)
%                04-April-2023 (TL, reorganized if/else)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%check if parantheses are used to select elements
if length(S) == 1 
    if strcmp(S.type,'()')
        %obtain sub-intervals from the interval object
        newObj = I;
        % only one index specified
        if length(S.subs)==1
            newObj.inf=I.inf(S.subs{1});
            newObj.sup=I.sup(S.subs{1});
        %two indices specified
        elseif length(S.subs)==2
            %Select column of obj
            if strcmp(S.subs{1},':')
                column=S.subs{2};
                newObj.inf=I.inf(:,column);
                newObj.sup=I.sup(:,column);
            %Select row of V    
            elseif strcmp(S.subs{2},':')
                row=S.subs{1};
                newObj.inf=I.inf(row,:);
                newObj.sup=I.sup(row,:);
            %Select single element of V    
            elseif isnumeric(S.subs{1})
                row=S.subs{1};
                column=S.subs{2};
                newObj.inf=I.inf(row,column);
                newObj.sup=I.sup(row,column);
            end
        end
        return
    elseif strcmp(S.type,'.') 
        if strcmp(S.subs,'inf')
            newObj = I.inf;
            return
        elseif strcmp(S.subs,'sup')
            newObj = I.sup;
            return
        end
    end   
end

% call build in subsref function as a default
newObj = builtin('subsref', I, S);

end

% ------------------------------ END OF CODE ------------------------------
