function obj = subsasgn(obj, S, value)
% subsasgn - Overloads the opertor that writes elements, e.g. I(1,2)=value,
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    obj = subsasgn(obj, S, value)
%
% Inputs:
%    obj - interval object 
%    S - contains information of the type and content of element selections
%    value - value to be written
%
% Outputs:
%    obj - interval object 
%
% Example: 
%    a=interval([-1 1], [1 2]);
%    a(1,2) = interval(-10,10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-June-2015 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%check if value is an interval
if ~isa(value,'interval')
    value = interval(value,value);
end

if ~isa(obj,'interval')
    obj = interval();
end

%check if parantheses are used to select elements
if strcmp(S.type,'()')
    % only one index specified
    if length(S.subs)==1
        obj.inf(S.subs{1}) = value.inf;
        obj.sup(S.subs{1}) = value.sup;
    %two indices specified
    elseif length(S.subs)==2
        %Select column of obj
        if strcmp(S.subs{1},':')
            column=S.subs{2};
            obj.inf(:,column) = value.inf;
            obj.sup(:,column) = value.sup;
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row=S.subs{1};
            obj.inf(row,:) = value.inf;
            obj.sup(row,:) = value.sup;
        %Select single element of V    
        elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
            row=S.subs{1};
            column=S.subs{2};
            obj.inf(row,column) = value.inf;
            obj.sup(row,column) = value.sup;
        end
    end
end

%------------- END OF CODE --------------