function pZ = deleteZeros(pZ)
% deleteZeros - (DEPRECATED -> compact)
%
% Syntax:
%    pZ = deleteZeros(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       20-April-2020 
% Last update:   ---
% Last revision: 30-July-2023 (MW, merged into compact)

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning("CORA:deprecated",'function','polyZonotope/deleteZeros','CORA v2024', ...
    'When updating the code, please replace every function call ''deleteZeros(pZ)'' with ''compact(pZ,''zeros'')''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')
pZ = compact_(pZ,'zeros',eps);

% ------------------------------ END OF CODE ------------------------------
