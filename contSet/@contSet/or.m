function S = or(varargin)
% or - computes an over-approximation for the union of two sets
%
% Description:
%    computes the set { s | s \in \mathcal{S}_1 \vee s \in \mathcal{S}_2 }
%
% Syntax:
%    S = S1 | S2
%    S = or(S1, S2)
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
