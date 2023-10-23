function S = and_(varargin)
% and_ - overloads '&' operator, computes the intersection of two sets
%    (internal use, see also contSet/and)
%
% Syntax:
%    S = and_(S1, S2)
%    S = and_(S1, S2, method)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%    method - (optional) char
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
