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

funcname = mfilename;
warning(sprintf(['The function ''' funcname ''' is deprecated (since CORA 2024) and has been replaced by ''compact''.\n' ...
    '         When updating the code, please rename every function call ''' funcname '(Z)'' -> ''compact(Z,''aligned'')''.\n' ...
    '         Note that the function ''' funcname ''' will be removed in a future release.']));
Z = compact_(Z,'aligned',eps);

% ------------------------------ END OF CODE ------------------------------
