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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2015 
% Last update:   13-October-2024 (MW, simplify, extend to nD intervals)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if value is an interval
if ~isa(val,'interval')
    val = interval(val,val);
end

if ~isa(I,'interval')
    I = interval.empty(dim(val));
end

% check if parentheses are used to select elements
if strcmp(S.type,'()')
    I.inf = builtin('subsasgn',I.inf,S,val.inf);
    I.sup = builtin('subsasgn',I.sup,S,val.sup);
end

% ------------------------------ END OF CODE ------------------------------
