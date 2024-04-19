function res = isHyperplane(S)
% isHyperplane - (DEPRECATED -> representsa)
%
% Syntax:
%    res = isHyperplane(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning('CORA:deprecated','function','isHyperplane','CORA v2024', ...
    'When updating the code, please replace every function call ''isHyperplane(S)'' with ''representsa(S,''hyperplane'')''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')
res = representsa(S,'hyperplane');

% ------------------------------ END OF CODE ------------------------------
