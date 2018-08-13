function [element] = subsref(V, S)
% subsref - Overloads the opertor that selects elements, e.g. V(1,2),
% where the element of the first row and second column is referred to.
%
% Syntax:  
%    [element] = subsref(V, S)
%
% Inputs:
%    V - vertices object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    V - elements of the vertices matrix
%
% Example: 
%    V=vertices(rand(2,6));
%    V(2,3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  31-January-2011
%               24-March-2015
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain vertices matrix from the vertices object
Vmatrix=V.V;

if ~isempty(Vmatrix)
    %check if parantheses are used to select elements
    if strcmp(S.type,'()')
        if length(S.subs)==2
            %Select column of V
            if strcmp(S.subs{1},':')
                column=S.subs{2};
                element=Vmatrix(:,column);
            %Select row of V    
            elseif strcmp(S.subs{2},':')
                row=S.subs{1};
                element=Vmatrix(row,:);
            %Select single element of V    
            elseif isnumeric(S.subs{1})&isnumeric(S.subs{1})
                row=S.subs{1};
                column=S.subs{2};
                element=Vmatrix(row,column);
            end
        %no selection if elements not proper specified  
        else
            element=[];
        end
    %no selection if parantheses are not used    
    else
        element=[];
    end

%     %write elements to V
%     V.V=element;
end

%------------- END OF CODE --------------