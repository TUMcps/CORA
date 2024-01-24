function res = in(S1,S2,varargin)
% in - (DEPRECATED -> contains)
%
% Syntax:
%    res = contains(S1,S2,tol)
%    res = in(S1,S2,tol) % deprecated
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contains

% Authors:       Mark Wetzlinger
% Written:       25-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

warning(sprintf(['The function ''in'' is deprecated (since CORA 2022) and has been replaced by ''contains''.\n' ...
    '         The main reason is that the syntax in(S1,S2) can also be written as S1.in(S2).\n' ...
    '         This is, however, the opposite of the actual computation, which checks whether S2 is a subset of S1.\n' ...
    '         The function ''contains'' eliminates this potential confusion, as S1.contains(S2) is semantically correct.\n' ...
    '         When updating the code, please rename every function call ''in(S1,S2)'' -> ''contains(S1,S2)''.\n' ...
    '         Note that the function ''in'' will be removed in a future release.']));
res = contains(S1,S2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
