function res = isInterval(pZ)
% isInterval - Checks if a polynomial zonotope represents an interval
%
% Syntax:  
%    res = isInterval(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([-0.5;0],[1.5 0;0 -2],[],[0 1;1 0]);
%    pZ2 = polyZonotope([-0.5;0],[-0.5 -0.5;0.5 -2],[],[1 1;0 1]);
%   
%    isZonotope(pZ1)
%    isZonotope(pZ2)
%
%    figure; hold on;
%    plot(pZ1,[1,2],'FaceColor','b','Splits',10);
% 
%    figure; hold on;
%    plot(pZ2,[1,2],'FaceColor','r','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, isPolytope, isZonotope

% Author:       Niklas Kochdumper
% Written:      14-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

% check if polynomial zonotope is a zonotope
if ~isZonotope(pZ)
    return;
end

% convert to zonotope and check is an interval
res = isInterval(zonotope(pZ));

%------------- END OF CODE --------------