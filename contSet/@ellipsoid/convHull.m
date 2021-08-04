function E = convHull(obj,S)
% convHull - computes an overapproximation of the convex hull of an
% ellipsoid and contSet
%
% Syntax:  
%    E = or(obj,S)
%
% Inputs:
%    obj                - ellipsoid object
%      S                - set (or cell array)
%
% Outputs:
%    E - resulting ellipsoid
%
% Example: 
%    E1 = ellipsoid.generateRandom(2,false);
%    E2 = ellipsoid.generateRandom(2,false); 
%    E = convHull(E1,E2);
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/or

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  19-March-2021 (complete rewrite)
% Last revision:---

%------------- BEGIN CODE --------------
% simply call "or"
E = or(obj,S,'o');
%------------- END OF CODE --------------