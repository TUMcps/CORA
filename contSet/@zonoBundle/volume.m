function [vol] = volume(Zbundle)
% volume - Computes the volume of a zonotope bundle
%
% Syntax:  
%    [vol] = volume(Zbundle)
%
% Inputs:
%    Zbundle - zonotope bundle
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain polytope of zonotope bundle
P = polytope(Zbundle);

%compute volume
vol = volume(P);


%------------- END OF CODE --------------
