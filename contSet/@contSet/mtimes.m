function S = mtimes(varargin)
% mtimes - overloaded '*' operator for the multiplication of a matrix
%    with a set
%
% Description:
%    computes the set { M s | s \in \mathcal{S} }, 
%    \mathcal{S} \in \mathbb{R}^n, M \in \mathbb{R}^{w \times n}
%
% Syntax:
%    S = M * S
%    S = mtimes(M, S)
%
% Inputs:
%    M - numeric
%    S - contSet object
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
