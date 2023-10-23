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

funcname = mfilename;
warning(sprintf(['The function ''' funcname ''' is deprecated (since CORA 2024) and has been replaced by ''compact''.\n' ...
    '         When updating the code, please rename every function call ''' funcname '(Z)'' -> ''compact(Z,''zeros'')''.\n' ...
    '         Note that the function ''' funcname ''' will be removed in a future release.']));
Z = compact_(Z,'zeros',eps);

% ------------------------------ END OF CODE ------------------------------
