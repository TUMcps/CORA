function res = isequal(O,S,varargin)
% isequal - checks if an empty set is equal to another set or point
%
% Syntax:
%    res = isequal(O,S)
%    res = isequal(O,S,tol)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    O1 = emptySet(2);
%    O2 = emptySet(3);
%    res1 = isequal(O1,O1);
%    res2 = isequal(O1,O2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% default values
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{O,'att',{'emptySet','numeric'}};
                {S,'att',{'contSet','numeric'}};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% ensure that numeric is second input argument
[O,S] = reorderNumeric(O,S);

% ambient dimensions must match
if ~equalDimCheck(O,S,true)
    res = false;
    return
end

if isa(S,'emptySet')
    % note: tolerance has no effect, only for overloading purposes
    res = O.dimension == S.dimension;
    return
end

% other contSet classes
if isa(S,'contSet')
    res = dim(S) == O.dimension && representsa_(S,'emptySet',eps);
    return
end

% vector
if isnumeric(S)
    res = isempty(S) && numel(S) == O.dimension;
    return
end

throw(CORAerror('CORA:noops',O,S));

% ------------------------------ END OF CODE ------------------------------
