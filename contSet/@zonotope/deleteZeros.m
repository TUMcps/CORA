function Z = deleteZeros(Z)
% deleteZeros - (DEPRECATED -> compact)
%
% Syntax:
%    Z = deleteZeros(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/compact_

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       15-January-2009
% Last update:   27-August-2019
% Last revision: 29-July-2023 (MW, merged into compact)

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated",'function','zonotope/deleteZeros','CORA v2024', ...
    'When updating the code, please replace every function call ''deleteZeros(Z)'' with ''compact(Z,''zeros'')''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')
Z = compact_(Z,'zeros',eps);

% ------------------------------ END OF CODE ------------------------------
