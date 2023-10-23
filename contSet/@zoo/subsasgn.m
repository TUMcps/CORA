function obj = subsasgn(obj, S, value)
% subsasgn - Overloads the operator that writes elements, e.g. T(1,2)=value,
% where the element of the first row and second column is referred to.
%
% Syntax:
%    obj = subsasgn(obj, S, value)
%
% Inputs:
%    obj - taylm object 
%    S - contains information of the type and content of element selections
%    value - value to be written
%
% Outputs:
%    obj - taylm object 
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Authors:       Niklas Kochdumper
% Written:       10-April-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{obj,'att','taylm'};
                {S,'att',{'numeric','struct'}};
                {value,'att','numeric','nonnan'}});

% loop over all combined calls to subsref
for j = 1:length(S)
    
    % override subsref for type obj(i,j)
    if strcmp(S(j).type,'()')
        
        % loop over all objects stored in the zoo
        for i = 1:length(obj.method)
            obj.objects{i} = subsasgn(obj.objects{i},S(j),value.objects{i});
        end
        
    % for all other types, call the build in function
    else
        obj = builtin('subsasgn', obj, S(j), value);
    end
end

% ------------------------------ END OF CODE ------------------------------
