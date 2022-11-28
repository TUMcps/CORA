function P = mptPolytope(hs)
% mptPolytope - Converts a halfspace to a polytope
%
% Syntax:  
%    P = mptPolytope(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    P - polytope object
%
% Example:
%    hs = halfspace([1 1],2);
%    P = mptPolytope(hs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/mptPolytope

% Author:       Victor Charlent
% Written:      28-June-2016
% Last update:  17-March-2017, Matthias Althoff
% Last revision:---

%------------- BEGIN CODE --------------

P = mptPolytope(hs.c',hs.d);

%------------- END OF CODE --------------
