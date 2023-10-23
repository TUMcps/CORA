function S = minkDiff(varargin)
% minkDiff - Minkowski difference
%
% Description:
%    computes the set { s | s \oplus \mathcal{S}_2 \subseteq \mathcal{S}_1 }
%
% Syntax:
%    S = minkDiff(S1,S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object, numerical vector
%    method - (optional) method
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       02-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:noops',varargin{:}));

% ------------------------------ END OF CODE ------------------------------
