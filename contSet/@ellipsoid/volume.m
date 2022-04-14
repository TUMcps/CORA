function val = volume(E)
% volume - Computes the volume of the ellipsoid E according to Sec. 2 in
%          [1]
%
% Syntax:  
%    val = volume(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    val - volume
%
% Example: 
%    E=ellipsoid([1,0;0,3],[1;-1]);
%    val = volume(E);
%
% References:
%    [1] A. Moshtagh. "Minimum volume enclosing ellipsoid", 2005
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      28-August-2019
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------
if isempty(E)
    val = 0;
    return;
end
n = length(E.Q);
% use volume for n-ball
Vball = pi^(n/2)/gamma(n/2+1);
val = Vball*sqrt(det(E.Q));
%------------- END OF CODE --------------