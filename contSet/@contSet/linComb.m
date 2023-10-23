function S = linComb(varargin)
% linComb - computes the linear combination of two sets
%
% Description:
%    computes the set { \lambda s_1 + (1-\lambda) s_2 | s_1 \in \mathcal{S}_1, s_2 \in \mathcal{S}_2, \lambda \in [0,1] }
%
% Syntax:
%    S = linComb(S1, S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
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
