function S = plus(varargin)
% plus - overloaded '+' operator for the Minkowski addition of two
%    sets or a set with a vector
%
% Description:
%    computes the set { s_1 + s_2 | s_1 \in \mathcal{S}_1, s_2 \in \mathcal{S}_2 }
%
% Syntax:
%    S = S1 + S2
%    S = plus(S1,S2)
%
% Inputs:
%    S1 - contSet object or numeric
%    S2 - contSet object or numeric
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
