function res = isIntersecting_(varargin)
% isIntersecting_ - checks if two sets intersect
%    (internal use, see also contSet/isIntersecting)
%
% Syntax:
%    S = isIntersecting_(S1, S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet objects
%    type - type of check 
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
