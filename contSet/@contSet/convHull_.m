function S = convHull_(varargin)
% convHull_ - computes an enclosure for the convex hull of a set and
%    another set or a point
%
% Description:
%    computes the set { \lambda s_1 + (1-\lambda) s_2 | s_1,s_2 \in \mathcal{S}_1 \cup \mathcal{S}_2, \lambda \in [0,1] }
%
% Syntax:
%    S = convHull_(S1, S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object or numeric
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: 30-September-2024 (MW, renamed convHull_)

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
