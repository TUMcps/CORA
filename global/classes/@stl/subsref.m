function res = subsref(obj,S)
% subsref - Overloads the operator that selects elements, e.g. obj(2)
%
% Syntax:  
%    res = subsref(obj,S)
%
% Inputs:
%    obj - stl object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    res - selected element represented as an stl object 
%
% Example: 
%    x = stl('x',2)
%    x(2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % check if parantheses are used to select elements
    if length(S) == 1 && strcmp(S.type,'()') 
        if length(S.subs) == 1     
            res = obj;
            res.var = obj.variables{S.subs{1}};
        else
            throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
        end
    else
        res = builtin('subsref',obj,S);
    end
end

%------------- END OF CODE --------------