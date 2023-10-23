function S = quadMap(varargin)
% quadMap - computes the quadratic map of a set
%
% Syntax:
%    S = quadMap(S1, S2, Q)
%
% Description:
%    computes the set { x | x_{(i)} = s_1^T Q_i s_2, ~s_1 \in \mathcal{S}_1,~s_2 \in \mathcal{S}_2, i = 1 \dots w \}, 
%    Q_i \in \R^{n \times n}
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%    Q - quadratic coefficients as a cell of matrices
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
