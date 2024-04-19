function res = isParallelotope(S)
% isParallelotope - (DEPRECATED -> representsa)
%
% Syntax:
%    res = isParallelotope(S)
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

CORAwarning('CORA:deprecated','function','contSet/isParallelotope','CORA v2024', ...
    'When updating the code, please replace every function call ''isParallelotope(S)'' with ''representsa(S,''parallelotope'')''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')
res = representsa(S,'parallelotope');

% ------------------------------ END OF CODE ------------------------------
