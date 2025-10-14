function newObj = subsref(matZ, S)
% subsref - Overloads the opertor that selects elements, e.g. I(1,2),
%    where the element of the first row and second column is referred to.
%
% Syntax:
%    newObj = subsref(matZ,S)
%
% Inputs:
%    matZ - matZonotope object 
%    S - contains information of the type and content of element selections  
%
% Outputs:
%    newObj - matZonotope object 
%
% Example: 
%    C = [ -0.537 -1.980 0.849 ; 0.329 -1.867 0.405 ; 1.054 -1.832 -0.702 ];
%    G{1} = [ 1.499 -1.019 -0.601 ; 0.138 -1.385 -1.172 ; -1.587 0.955 -0.577 ];
%    matZ = matZonotope(C,G);
%    matZ(1,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper, Tobias Ladner
% Written:       09-November-2018 
% Last update:   12-November-2018 (NK, default to build in for other cases)
%                25-April-2024 (TL, faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%check if parantheses are used to select elements
if length(S) == 1 && strcmp(S.type,'()')
    
    %obtain sub-intervals from the interval object
    newObj = matZ;

    [n,m,h] = size(matZ.G);
    
    % only one index specified
    if length(S.subs)==1
        newObj.C = matZ.C(S.subs{1});
        G = reshape(matZ.G, n*m,h);
        newObj.G = G(S.subs{1},:);
    %two indices specified
    elseif length(S.subs)==2
        %Select column of obj
        if strcmp(S.subs{1},':')
            column = S.subs{2};
            newObj.C = matZ.C(:,column);
            newObj.G = matZ.G(:,column,:);
        %Select row of V    
        elseif strcmp(S.subs{2},':')
            row = S.subs{1};
            newObj.C = matZ.C(row,:);
            newObj.G = matZ.G(row,:,:);
        %Select single element of V    
        elseif isnumeric(S.subs{1}) && isnumeric(S.subs{2})
            row = S.subs{1};
            column = S.subs{2};
            newObj.C = matZ.C(row,column);
            newObj.G = matZ.G(row,column);
        end
    end
else
    % call build in subsref function as a default
    newObj = builtin('subsref', matZ, S);
end

% ------------------------------ END OF CODE ------------------------------
