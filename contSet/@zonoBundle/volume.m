function vol = volume(zB)
% volume - Computes the volume of a zonotope bundle
%
% Syntax:  
%    vol = volume(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    vol - volume
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-February-2011 
% Last update:  18-August-2022 (MW, include standardized preprocessing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[vol,vars] = pre_volume('zonoBundle',zB);

% check premature exit
if vol
    % if result has been found, it is stored in the first entry of var
    vol = vars{1}; return
end

%obtain polytope of zonotope bundle
P = mptPolytope(zB);

%compute volume
vol = volume(P);

%------------- END OF CODE --------------
