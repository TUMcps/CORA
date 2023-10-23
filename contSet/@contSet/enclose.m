function S = enclose(varargin)
% enclose - encloses a set and its affine transformation
%
% Description:
%    Computes the set
%    { \lambda x1 + (1 - \lambda) * x2 | x1 \in S1, x2 \in S2, \lambda \in [0,1] }
%    where S2 = M*S1 + Splus
%
% Syntax:
%    S = enclose(S1, S2)
%    S = enclose(S1, M, Splus)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%    Splus - contSet object
%    M - matrix for the linear transformation
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
