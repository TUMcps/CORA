function newObj = subsref(I,S)
% subsref - Overloads the operator that selects elements, e.g., I(1,2),
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
%    I = interval([-1; 1], [1; 2]);
%    I(2,1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       19-June-2015 
% Last update:   22-June-2015 
%                12-November-2018 (NK, default to build in for other cases)
%                04-April-2023 (TL, reorganized if/else)
%                13-October-2024 (MW, simplify using builtin-subsref)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%check if parentheses are used to select elements
if length(S) == 1
    if strcmp(S.type,'()')
        % obtain sub-intervals from the interval object
        newObj = I;
        newObj.inf = builtin('subsref', newObj.inf, S);
        newObj.sup = builtin('subsref', newObj.sup, S);
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

% call built-in subsref function as a default
newObj = builtin('subsref', I, S);

% ------------------------------ END OF CODE ------------------------------
