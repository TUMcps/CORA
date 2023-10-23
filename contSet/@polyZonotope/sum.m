function pZsum = sum(pZ,pZlist)
% sum - computes the sum of multiple polynomial zonotopes
%
% Syntax:
%    pZsum = sum(pZ,pZlist)
%
% Inputs:
%    pZ - polyZonotope object (first part of the sum)
%    pZlist - list of polyZonotope objects (cell-array)
%
% Outputs:
%    pZsum - polyZonotope object representing the set of the sum
%
% Example: 
%    pZ1 = polyZonotope([1;2],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[2 0 1;0 2 1],[],[0 0 0;0 0 0;1 0 3;0 1 1]);
%
%    pZsum = sum(pZ1,{pZ2,pZ1});
%
%    figure
%    plot(pZ1,[1,2],'FaceColor','r');
%    figure
%    plot(pZ2,[1,2],'FaceColor','b');
%    figure
%    plot(pZsum,[1,2],'FaceColor','g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Niklas Kochdumper
% Written:       17-August-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize variables
EList = cell(length(pZlist)+1,1);
G = pZ.G;
GI = pZ.GI;
c = pZ.c;

% bring exponent matrices to a common representation
EList{1} = pZ.E;
id = pZ.id;

for i = 1:length(pZlist)
    [id,~,EList{i+1}] = mergeExpMatrix(id,pZlist{i}.id, ...
                                    EList{i},pZlist{i}.E);
end

% concatenate exponent and generator matrices and sum up center
M = length(id);
E = [pZ.E;zeros(M-size(pZ.E,1),size(pZ.E,2))];

for i = 1:length(pZlist)
   G = [G,pZlist{i}.G];
   GI = [GI,pZlist{i}.GI];
   eTemp = EList{i+1};
   E = [E,[eTemp;zeros(M-size(eTemp,1),size(eTemp,2))]];
   c = c + pZlist{1}.c;
end

% add up all generators that belong to identical exponents
[ExpNew,Gnew] = removeRedundantExponents(E,G);

% construct the resulting polynomial zonotope
pZsum = polyZonotope(c,Gnew,GI,ExpNew,id);

% ------------------------------ END OF CODE ------------------------------
