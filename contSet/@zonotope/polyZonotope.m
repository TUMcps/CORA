function pZ = polyZonotope(Z)
% polyZonotope - converts a zonotope object to a polyZonotope object
%
% Syntax:
%    pZ = polyZonotope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    Z = zonotope([1 2 0 -1;3 1 2 2]);
%    pZ = polyZonotope(Z);
%    
%    figure; xlim([-3,5]); ylim([-3,9]);
%    plot(Z,[1,2],'FaceColor','b');
%    plot(pZ,[1,2],'r--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/polyZonotope, taylm/polyZonotope

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = Z.c;
G = Z.G;
E = eye(size(G,2));

pZ = polyZonotope(c,G,[],E);
    
% ------------------------------ END OF CODE ------------------------------
