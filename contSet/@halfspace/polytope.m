function P = polytope(hs)
% polytope - Converts a halfspace to a polytope
%
% Syntax:
%    P = polytope(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    P - polytope object
%
% Example:
%    hs = halfspace([1 1],2);
%    P = polytope(hs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/polytope

% Authors:       Victor Charlent
% Written:       28-June-2016
% Last update:   17-March-2017 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

P = polytope(hs.c',hs.d);

% set properties
P.bounded.val = false;

% ------------------------------ END OF CODE ------------------------------
