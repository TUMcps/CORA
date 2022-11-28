function t = taylm(Z)
% taylm - enclose a zonotope object with a Taylor model
%
% Syntax:  
%    t = taylm(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    t - taylm object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% create taylor models for factors
m = size(Z.Z,2)-1;
dom = interval(-ones(m,1),ones(m,1));
t = taylm(dom);

% create taylor model for the zonotope
t = Z.Z(:,1) + Z.Z(:,2:end)*t;

%------------- END OF CODE --------------