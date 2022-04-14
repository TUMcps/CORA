function newObj = subsref(obj, S)
% subsref - Overloads the opertor that selects elements, e.g. I(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    newObj = subsref(obj, S)
%
% Inputs:
%    obj - matrix zonotope object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    newObj - matrix zonotope object 
%
% Example: 
%    C = rand(3);
%    G{1} = rand(3);
%    matZ = matZonotope(C, G);
%    matZ(1,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      09-November-2018 
% Last update:  12-November-2018 (NK: default to build in for other cases)
% Last revision:---

%------------- BEGIN CODE --------------

%check if parantheses are used to select elements
if length(S) == 1 && strcmp(S.type,'()')
    
    %obtain sub-intervals from the interval object
    newObj = obj;
    
    % only one index specified
    if length(S.subs)==1
        newObj.center = obj.center(S.subs{1});
        for iGen = 1:newObj.gens
            newObj.generator{iGen} = obj.generator{iGen}(S.subs{1});
        end
    %two indices specified
    elseif length(S.subs)==2
        %Select column of obj
        if strcmp(S.subs{1},':')
            column = S.subs{2};
            newObj.center = obj.center(:,column);
            for iGen = 1:newObj.gens
                newObj.generator{iGen} = obj.generator{iGen}(:,column);
            end
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row = S.subs{1};
            newObj.center = obj.center(row,:);
            for iGen = 1:newObj.gens
                newObj.generator{iGen} = obj.generator{iGen}(row,:);
            end
        %Select single element of V    
        elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
            row = S.subs{1};
            column = S.subs{2};
            newObj.center = obj.center(row,column);
            for iGen = 1:newObj.gens
                newObj.generator{iGen} = obj.generator{iGen}(row,column);
            end
        end
    end
else
    % call build in subsref function as a default
    newObj = builtin('subsref', obj, S);
end

%------------- END OF CODE --------------