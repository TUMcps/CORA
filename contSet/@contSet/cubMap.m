function S = cubMap(varargin)
% cubMap - computes the quadratic map of a set
%
% Description:
%    Calculates the following set:
%    { x | x_{(i)} = \sum_{j=1}^n s_{(j)} ~(s^T T_{i,j} s), s \in \mathcal{S}, i = 1 \dots w \},
%    T_{i,j} \in \R^{n \times n}
%
% Syntax:
%    res = cubMap(S1,T)
%    res = cubMap(S1,T,ind)
%    res = cubMap(S1,S2,S3,T)
%    res = cubMap(S1,S2,S3,T,ind)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%    S3 - contSet object
%    T - third-order tensor
%    ind - cell-array containing the non-zero indices of the tensor
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
