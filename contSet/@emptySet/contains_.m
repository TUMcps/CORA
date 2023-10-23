function res = contains_(O,S,type,tol,varargin)
% contains_ - determines if an empty set contains a set or a point
%
% Syntax:
%    res = contains_(O,S)
%    res = contains_(O,S,type)
%    res = contains_(O,S,type,tol)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    type - 'exact' or 'approx'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    O = emptySet(2);
%    p = [1;1];
%    res = contains(O,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, rename contains_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimensions are already checked...
res = false;
if isa(S,'emptySet')
    % empty set contains the empty set
    res = true;

elseif isa(S,'contSet') && representsa_(S,'emptySet',eps)
    % empty set contains contSet objects if they also represent the empty
    % set
    res = true;

elseif isnumeric(S) && isempty(S)
    % empty set contains empty vectors
    res = true;

end

% ------------------------------ END OF CODE ------------------------------
