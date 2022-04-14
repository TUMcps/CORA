function pZ = polyZonotope(obj)
% polyZonotope - Convert zonotope bundle to a polynomial zonotope
%
% Syntax:  
%    pZ = polyZonotope(obj)
%
% Inputs:
%    obj - zonoBundle object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    pZ = polyZonotope(zB);
%
%    figure
%    plot(zB,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-1,4]);
%    ylim([-4,3]);
% 
%    figure
%    plot(pZ,[1,2],'b','Filled',true,'EdgeColor','none','Splits',12);
%    xlim([-1,4]);
%    ylim([-4,3]);
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