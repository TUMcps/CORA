function res = contains_(cPZ,S,type,varargin)
% contains_ - determines if a constrained polynomial zonotope contains a set
%    or a point
%
% Syntax:
%    res = contains_(cPZ,S)
%    res = contains_(cPZ,S,type)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object or single point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    E = [1 0 3; 0 1 1; 0 0 0];
%    GI = [0.5;0]; 
%    A = [1 1 -1.5];
%    b = 0.5;
%    EC = [1 0 0;0 1 0;0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
% 
%    c = [0;0];
%    G = [1 -2 1; 2 3 1];
%    E = [1 0 2;0 1 1; 0 0 0];
%    A = [1 -1 0.5];
%    b = -0.5;
%    EC = [2 0 0;0 1 0; 0 0 1];
%    cPZ1 = conPolyZono(c,0.3*G,E,A,b,EC);
%    cPZ2 = conPolyZono(c,0.5*G,E,A,b,EC);
% 
%    p1 = [1;1];
%    p2 = [-1;3];
% 
%    contains(cPZ,p1,'approx')
%    contains(cPZ,p2,'approx')
%    contains(cPZ,cPZ1,'approx')
%    contains(cPZ,cPZ2,'approx')
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Splits',12);
%    plot(p1(1),p1(2),'.g','MarkerSize',20);
%    plot(p2(1),p2(2),'.r','MarkerSize',20);
%    plot(cPZ1,[1,2],'g','Splits',12);
%    plot(cPZ2,[1,2],'r','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains_, conZonotope/contains_

% Authors:       Niklas Kochdumper
% Written:       05-February-2020 
% Last update:   25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% check user inputs 
if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',cPZ,S));
end

% transform to equivalent higher-dimensional polynomial zonotope
m = size(cPZ.A,1);
c = [cPZ.c; -cPZ.b];
G = blkdiag(cPZ.G,cPZ.A);
E = [cPZ.E,cPZ.EC];

GI = cPZ.GI;
if ~isempty(GI)
    GI = [GI; zeros(m,size(GI,2))];
end

pZ = polyZonotope(c,G,GI,E,cPZ.id);

% increase dimension of the second set to match dim. of first set
if ~isempty(cPZ.A)
    if isnumeric(S)
        S = [S;ones(m,1)];
    else
        temp = sqrt(max(sum(cPZ.A.^2,1)))/100 * ones(m,1);
        S = cartProd_(S,interval(-temp,temp),'exact');
    end
end

% check containment using the function for higher-dimensional zonotopes
res = contains_(pZ,S,'approx');

% ------------------------------ END OF CODE ------------------------------
