function vol = volume_(E,varargin)
% volume_ - computes the volume of an ellipsoid acc. to Sec. 2 in [1]
%
% Syntax:
%    vol = volume_(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    vol - volume
%
% Example: 
%    E = ellipsoid([1,0;0,3],[1;-1]);
%    vol = volume(E);
%
% References:
%    [1] A. Moshtagh. "Minimum volume enclosing ellipsoid", 2005
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Victor Gassmann
% Written:       28-August-2019
% Last update:   04-July-2022 (VG, allow class array input)
%                18-August-2022 (MW, include standardized preprocessing)
%                05-October-2024 (MW, remove class array)
% Last revision: 27-March-2023 (MW, rename volume_)

% ------------------------------ BEGIN CODE -------------------------------

% empty ellipsoid
if representsa_(E,'emptySet',eps)
    vol = 0;
    return
end

% dimension
n = dim(E);

% use volume for n-ball
Vball = pi^(n/2)/gamma(n/2+1);
vol = Vball*sqrt(det(E.Q));

% ------------------------------ END OF CODE ------------------------------
