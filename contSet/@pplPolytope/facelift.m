function [P] = facelift(P,Z)
% facelift - face lifting such that the enlarged polytope encloses a
% Minkowski addition with the zonotope Z
%
% Syntax:  
%    [P] = facelift(P,Z)
%
% Inputs:
%    P - pplPolytope object
%    Z - zonotope object
%
% Outputs:
%    P - pplPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      02-March-2011
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%enlarge for each halfspace
for iHalfspace = 1:length(P.d)
    d_IH = interval(P.C(iHalfspace,:)*Z);
    d_add = supremum(d_IH);
    P.d(iHalfspace) = P.d(iHalfspace) + d_add;
end


%------------- END OF CODE --------------