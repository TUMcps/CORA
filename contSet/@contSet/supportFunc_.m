function S = supportFunc_(varargin)
% supportFunc_ - evaluates the support function of a set along a given direction
%    (internal use, see also contSet/supportFunc)
%
% Syntax:
%    val = supportFunc_(S, dir)
%    val = supportFunc_(S, dir, type)
%    val = supportFunc_(S, dir, type, method, varargin)
%
% Inputs:
%    S - contSet object
%    dir - column vector specifying the direction 
%          for which the bounds are calculated
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    method - (optional) method that is used to calculate the bounds
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
