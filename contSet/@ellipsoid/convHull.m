function E = convHull(E,S)
% convHull - computes an overapproximation of the convex hull of an
%    ellipsoid and another set representation
%
% Syntax:  
%    E = or(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet class object (or class array)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1 = ellipsoid.generateRandom('Dimension',2);
%    E2 = ellipsoid.generateRandom('Dimension',2);
%    E = convHull(E1,E2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/or

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  19-March-2021 (complete rewrite)
%               04-July-2022 (input checks)
% Last revision:---

%------------- BEGIN CODE --------------

% check inputs
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}});

% convex hull with empty set
if isempty(S)
    return;
end

% check dimensions
equalDimCheck(E,S);

% simply call "or"
E = or(E,S,'outer');

%------------- END OF CODE --------------