function cPZ = exactPlus(cPZ,cPZ2)
% exactPlus - exact plus of two constrained polynomial zonotopes
%
% Syntax:
%    cPZ = exactPlus(cPZ,cPZ2)
%
% Inputs:
%    cPZ - conPolyZono object
%    cPZ2 - conPolyZono object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    E = [1 0 3;0 1 1];
%    A = [1 -1];
%    b = 0;
%    EC = [2 0; 0 1];
%    cPZ1 = conPolyZono(c,G,E,A,b,EC);
%
%    M = 0.1 * [3 1;2 4];
%    cPZ2 = M * cPZ1;
%
%    res = exactPlus(cPZ1,cPZ2);
%    res_ = cPZ1 + cPZ2;
%
%    figure; hold on
%    h1 = plot(res_,[1,2],'FaceColor','b','Splits',10);
%    h2 = plot(res,[1,2],'FaceColor','r','Splits',10);
%    legend([h1;h2],'Minkowski sum','exact plus');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, polyZonotope/exactPlus

% Authors:       Niklas Kochdumper
% Written:       14-August-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if constraints are identical
[id,EC1,EC2] = mergeExpMatrix(cPZ.id,cPZ2.id,cPZ.EC,cPZ2.EC);

if any(size(cPZ.A)~=size(cPZ2.A)) || any(any(abs(cPZ.A-cPZ2.A) > eps)) || ...
   any(size(cPZ.A)~=size(cPZ2.A)) || any(any(abs(cPZ.b-cPZ2.b) > eps)) || ...``
   any(size(cPZ.A)~=size(cPZ2.A)) || any(any(abs(EC1-EC2) > eps))
  throw(CORAerror('CORA:specialError',...
      'Operation only defined for conPolyZonotopes with identical constraints!'));
end

% call exactPlus for polynomial zonotopes
S1 = polyZonotope(cPZ.c,cPZ.G,cPZ.GI, cPZ.E,cPZ.id);
S2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.GI, cPZ2.E,cPZ2.id);

pZ = exactPlus(S1,S2);

% construct resulting constrained polynomial zonotope
[id,E,EC] = mergeExpMatrix(pZ.id,id,pZ.E,EC1);

cPZ = conPolyZono(pZ.c,pZ.G,E,cPZ.A,cPZ.b,EC,pZ.GI,id);

% ------------------------------ END OF CODE ------------------------------
