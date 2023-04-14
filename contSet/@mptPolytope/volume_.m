function vol = volume_(P,varargin)
% volume_ - Computes the volume of a mptPolytope
%
% Syntax:  
%    vol = volume_(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    vol - volume
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-February-2011 
% Last update:  18-August-2022 (MW, include standardized preprocessing)
% Last revision:27-March-2023 (MW, rename volume_)

%------------- BEGIN CODE --------------

%call volume operation of mpt toolbox
vol = volume(P.P);

%------------- END OF CODE --------------