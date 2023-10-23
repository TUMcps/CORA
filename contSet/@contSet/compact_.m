function S = compact_(varargin)
% compact_ - removes redundancies in the representation of a set, the
%    resulting set is equal to the original set, but minimal in its
%    representation
%    (internal use, see also contSet/compact)
%
% Syntax:
%    S = compact_(S)
%    S = compact_(S, method, tol)
%
% Inputs:
%    S - contSet object
%    method - (optional) method
%    tol - (optional) tolerance
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
