function pZ = polyZonotope(obj)
% polyZonotope - Convert constrained zonotope to a polynomial zonotope
%
% Syntax:  
%    pZ = polyZonotope(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    pZ = polyZonotope(cZ);
%
%    figure
%    plot(cZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-4,5]);
%    ylim([-2,3]);
% 
%    figure
%    plot(pZ,[1,2],'b','Filled',true,'EdgeColor','none');
%    xlim([-4,5]);
%    ylim([-2,3]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, mptPolytope/polyZonotope

% Author:       Niklas Kochdumper
% Written:      26-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    poly = mptPolytope(obj);
    pZ = polyZonotope(poly);

%------------- END OF CODE --------------