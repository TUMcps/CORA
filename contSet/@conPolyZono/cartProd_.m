function cPZ = cartProd_(cPZ,S,varargin)
% cartProd_ - Returns the Cartesian product of a constrained polynomial
%    zonotope and another set
%
% Syntax:
%    cPZ = cartProd_(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%
% Outputs:
%    cPZ - conPolyZono object representing the Cartesian product
%
% Example: 
%    c = [0;0];
%    G = [2 2; 2 -1];
%    E = [1 0; 0 1; 0 0];
%    A = [1 1 -1];
%    b = 2;
%    EC = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    I = interval(1,2);
%
%    res = cartProd(cPZ,I);
%
%    figure; hold on; box on; grid on;
%    plot(res,[1,2,3],'FaceColor','b','Splits',10);
%    plot(cPZ,[1,2],'FaceColor','r','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, polyZonotope/cartProd_, zonotope/cartProd_

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename cartProd_)

% ------------------------------ BEGIN CODE -------------------------------

% convert other set representations to polyZonotopes (first set)
if ~isa(cPZ,'conPolyZono')
    if isa(cPZ,'zonotope') || isa(cPZ,'interval') || isnumeric(cPZ) 
    
        Z = zonotope(cPZ);
        cPZ = conPolyZono(center(Z),[],[],generators(Z));
        
    elseif isa(cPZ,'polytope') || isa(cPZ,'taylm') || ...
           isa(cPZ,'zonoBundle') || isa(cPZ,'conZonotope') || ...
           isa(cPZ,'ellipsoid') || isa(cPZ,'capsule') || ...
           isa(cPZ,'polyZonotope')
        
        cPZ = conPolyZono(cPZ);
   
    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',cPZ,S));
    end
end
    
if ~isa(S,'conPolyZono')
    if isa(S,'zonotope') || isa(S,'interval') || isnumeric(S)
        
        Z = zonotope(S);
        S = conPolyZono(center(Z),[],[],generators(Z));
        
    elseif isa(S,'polytope') || isa(S,'taylm') || ...
           isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
           isa(S,'ellipsoid') || isa(S,'capsule') || ...
           isa(S,'polyZonotope')
        
        S = conPolyZono(S);

    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',cPZ,S));
    end
end 
      
% compute Cartesian product for polynomial zonotopes
pZ1 = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
pZ2 = polyZonotope(S.c,S.G,S.GI,S.E,S.id);

pZ = cartProd_(pZ1,pZ2,'exact');

% update constraints
cPZ = updateConstraints(conPolyZono(pZ),cPZ,S);

% ------------------------------ END OF CODE ------------------------------
