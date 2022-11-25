function [vol] = volume(P)
% volume - Computes the volume of a mptPolytope
%
% Syntax:  
%    [vol] = volume(P)
%
% Inputs:
%    P - mptPolytope object
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

%call volume operation of mpt toolbox
vol = volume(P.P);


%------------- END OF CODE --------------