function res = isequal(O,S,varargin)
% isequal - checks if two emptySet objects are equal
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

if isa(S,'emptySet')
    % note: tolerance has no effect, only for overloading purposes
    res = O.dimension == S.dimension;

elseif isnumeric(S)
    % vector
    res = isempty(S);

elseif isa(S,'contSet')
    % contSet objects...
    res = dim(S) == O.dimension && representsa_(S,'emptySet',eps);
    
else
    throw(CORAerror('CORA:noops',O,S));
end

% ------------------------------ END OF CODE ------------------------------
