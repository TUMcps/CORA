function res = cartProd_(varargin)
% cartProd_ - computes the Cartesian product of two sets
%    (internal use, see also contSet/cartProd)
%
% Syntax:
%    S = cartProd_(S1, S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%
% Outputs:
%    res - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
