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

funcname = mfilename;
warning(sprintf(['The function ''' funcname ''' is deprecated (since CORA 2024) and has been replaced by ''compact''.\n' ...
    '         When updating the code, please rename every function call ''' funcname '(pZ)'' -> ''compact(pZ,''zeros'')''.\n' ...
    '         Note that the function ''' funcname ''' will be removed in a future release.']));
pZ = compact_(pZ,'zeros',eps);

% ------------------------------ END OF CODE ------------------------------
