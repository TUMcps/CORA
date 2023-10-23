function val = norm_(varargin)
% norm_ - compute the norm of a set
%
% Syntax:
%    val = norm_(S,type)
%
% Inputs:
%    Z - contSet object
%    type - (optional) which kind of norm (default: 2)
%
% Outputs:
%    val - norm value
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/norm

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
