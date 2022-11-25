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
%    res - 1 if set is a interval, 0 if not
%
% Example: 
%    pZ1 = polyZonotope([-0.5;0],[1.5 0;0 -2],[],[0 1;1 0]);
%    pZ2 = polyZonotope([-0.5;0],[-0.5 -0.5;0.5 -2],[],[1 1;0 1]);
%   
%    isZonotope(pZ1)
%    isZonotope(pZ2)
%
%    figure
%    hold on
%    plot(pZ1,[1,2],'b','Filled',true,'EdgeColor','none','Splits',10);
%
%    figure
%    hold on
%    plot(pZ2,[1,2],'r','Filled',true,'EdgeColor','none','Splits',10);
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

    res = 0;

    % check if polynomial zonotope is a zonotope
    if ~isZonotope(pZ)
        return;
    end

    % convert to zonotope
    zono = zonotope(pZ);

    % check if zonotope is an interval
    res = isInterval(zono);
end

%------------- END OF CODE --------------