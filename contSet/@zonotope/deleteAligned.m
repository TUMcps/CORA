function Z = deleteAligned(Z)
% deleteAligned - (DEPRECATED -> compact)
%
% Syntax:
%    Z = deleteAligned(Z)
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
% See also: zonotope/reduce

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       15-January-2009
% Last update:   27-August-2019
% Last revision: 29-July-2023 (MW, merged to compact)

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated",'function','zonotope/deleteAligned','CORA v2024', ...
    'When updating the code, please replace every function call ''deleteAligned(Z)'' with ''compact(Z,''aligned'')''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')
Z = compact_(Z,'aligned',eps);

% ------------------------------ END OF CODE ------------------------------
