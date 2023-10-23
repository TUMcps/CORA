function pZ = polyZonotope(cZ)
% polyZonotope - Converts a constrained zonotope to a polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    pZ = polyZonotope(cZ);
%
%    figure; hold on; xlim([-4,5]); ylim([-2,3]);
%    plot(cZ,[1,2],'FaceColor','r');
%    plot(pZ,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, polytope/polyZonotope

% Authors:       Niklas Kochdumper
% Written:       26-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

poly = polytope(cZ);
pZ = polyZonotope(poly);

% ------------------------------ END OF CODE ------------------------------
