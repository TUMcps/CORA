function res = taylm(obj)
% taylm - enclose a zonotope object with a Taylor model
%
% Syntax:  
%    res = taylm(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%    res - taylm object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % create taylor models for factors
    m = size(obj.Z,2)-1;
    dom = interval(-ones(m,1),ones(m,1));
    t = taylm(dom);
    
    % create taylor model for the zonotope
    res = obj.Z(:,1) + obj.Z(:,2:end)*t;

%------------- END OF CODE --------------