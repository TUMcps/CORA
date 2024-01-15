function I = subsasgn(I,S,val)
% subsasgn - Overloads the operator that writes elements, e.g., I(1,2)=val,
%    where the element of the first row and second column is referred to.
%
% Syntax:
%    I = subsasgn(I,S,val)
%
% Inputs:
%    I - interval object 
%    S - contains information of the type and content of element selections
%    val - value to be inserted
%
% Outputs:
%    I - interval object 
%
% Example: 
%    I = interval([-1 1], [1 2]);
%    I(1,2) = interval(-10,10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       26-June-2015 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%check if value is an interval
if ~isa(val,'interval')
    val = interval(val,val);
end

if ~isa(I,'interval')
    I = interval.empty(dim(val));
end

%check if parantheses are used to select elements
if strcmp(S.type,'()')
    % only one index specified
    if length(S.subs)==1
        I.inf(S.subs{1}) = val.inf;
        I.sup(S.subs{1}) = val.sup;
    %two indices specified
    elseif length(S.subs)==2
        %Select column of obj
        if strcmp(S.subs{1},':')
            column=S.subs{2};
            I.inf(:,column) = val.inf;
            I.sup(:,column) = val.sup;
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row=S.subs{1};
            I.inf(row,:) = val.inf;
            I.sup(row,:) = val.sup;
        %Select single element of V    
        elseif isnumeric(S.subs{1})
            row=S.subs{1};
            column=S.subs{2};
            I.inf(row,column) = val.inf;
            I.sup(row,column) = val.sup;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
