function S = representsa_(varargin)
% representsa_ - checks if a set can also be represented by a different set,
%    e.g., a special case
%    (internal use, see also contSet/representsa)
%
% Syntax:
%    res = representsa_(S, type, tol)
%    [res,S] = representsa_(S, type, tol)
%
% Inputs:
%    S - contSet object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - (optional) tolerance
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
