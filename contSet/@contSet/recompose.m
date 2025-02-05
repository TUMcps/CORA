function S_out = recompose(S)
% recompose - overload for global/.../recompose function
%
% Syntax:
%    S_out = recompose(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    S - contSet object
%
% Example:
%    S = zonotope(0);
%    S_out = recompose(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: global/.../recompose, decompose, priv_reach_decomp

% Authors:       Mark Wetzlinger
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% nothing to recompose...
S_out = S.copy();

% ------------------------------ END OF CODE ------------------------------
