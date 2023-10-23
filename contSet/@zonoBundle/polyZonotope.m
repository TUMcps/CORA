function pZ = polyZonotope(zB)
% polyZonotope - Convert zonotope bundle to a polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    pZ = polyZonotope(zB);
%
%    figure; hold on; xlim([-1,4]); ylim([-4,3]);
%    plot(zB,[1,2],'FaceColor','r');
%    plot(pZ,[1,2],'FaceColor','b','Splits',12);
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

poly = polytope(zB);
pZ = polyZonotope(poly);

% ------------------------------ END OF CODE ------------------------------
