function pZ = polyZonotope(cPZ)
% polyZonotope - enclose a constrained polynomial zonotope with a
%    polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object 
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:  
%    c = [0;0];
%    G = [1 0 1 -1; 0 2 1 2];
%    E = [1 2 1 0; 0 0 1 2; 0 0 0 0];
%    A = [1 1 0.5];
%    b = 0.5;
%    EC = [0 1 0;1 0 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    pZ = polyZonotope(cPZ);
%   
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(pZ);
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope, interval, zonotope

% Authors:       Niklas Kochdumper
% Written:       04-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% remove all constraints
cPZ = reduceConstraints(cPZ,0);

% represent the resulting object as a polynomial zonotope
pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
    
% ------------------------------ END OF CODE ------------------------------
