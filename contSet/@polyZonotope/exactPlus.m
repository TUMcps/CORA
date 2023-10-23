function pZ = exactPlus(pZ1,pZ2)
% exactPlus - computes the addition of two sets while preserving the
%    dependencies between the two sets
%
% Syntax:
%    pZ = exactPlus(pZ1,pZ2)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    pZ2 = [1 2;-1 1]*pZ1;
%   
%    pZ = pZ1 + pZ2;
%    pZ_ = exactPlus(pZ1,pZ2);
%
%    figure
%    subplot(1,2,1);
%    plot(pZ,[1,2],'FaceColor','r','Splits',10);
%    title('Minkowski Sum');
%    subplot(1,2,2);
%    plot(pZ_,[1,2],'FaceColor','b');
%    title('Exact Addition');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, zonotope/plus

% Authors:       Niklas Kochdumper
% Written:       26-March-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% bring the exponent matrices to a common representation
[id,E1,E2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.E,pZ2.E);

% add up all generators that belong to identical exponents
[Enew,Gnew] = removeRedundantExponents([E1,E2],[pZ1.G,pZ2.G]);

% assemble the properties of the resulting polynomial zonotope
pZ = pZ1;
pZ.G = Gnew;
pZ.E = Enew;

pZ.c = pZ1.c + pZ2.c;
pZ.GI = [pZ1.GI,pZ2.GI];
pZ.id = id;

% ------------------------------ END OF CODE ------------------------------
