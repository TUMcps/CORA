function pZ = cartProd_(pZ,S,varargin)
% cartProd_ - Returns the Cartesian product a polynomial zonotope and
%    another set representation
%
% Syntax:
%    pZ = cartProd_(pZ,S)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope(2,[1 3 1],[],[1,2,3]);
%    Z = zonotope([1,3]);
%
%    pZcart = cartProd(pZ,Z);
%
%    figure; xlim([1 8]); ylim([-3 5]);
%    plot(pZcart,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename cartProd_)

% ------------------------------ BEGIN CODE -------------------------------

% convert other set representations to polyZonotopes (first set)
if ~isa(pZ,'polyZonotope')
    if isa(pZ,'zonotope') || isa(pZ,'interval')
        
        pZ = zonotope(pZ);
        pZ = polyZonotope(center(pZ),[],generators(pZ),[]);
        
    elseif isa(pZ,'conPolyZono')
        
        S = conPolyZono(S);
        pZ = cartProd_(pZ,S,'exact'); 
        return;
        
    elseif isa(pZ,'polytope') || isa(pZ,'zonoBundle') || ...
           isa(pZ,'conZonotope')

        pZ = polyZonotope(pZ);
        
    elseif isnumeric(pZ)
        pZ = polyZonotope(pZ,[],[],[]);
    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',pZ,S));
    end
end

% convert other set representations to polyZonotopes (second set)
if ~isa(S,'polyZonotope')
    if isa(S,'zonotope') || isa(S,'interval')
        
        S = zonotope(S);
        S = polyZonotope(center(S),[],generators(S),[]);
        
    elseif isa(S,'conPolyZono')
        
        pZ = conPolyZono(pZ);
        pZ = cartProd_(pZ,S,'exact'); 
        return;
        
    elseif isa(S,'polytope') || isa(S,'zonoBundle') || ...
           isa(S,'conZonotope')

        S = polyZonotope(S);
        
    elseif isnumeric(S)
        S = polyZonotope(S,[],[],[]);
    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',pZ,S));
    end
end

% get dimensions
n1 = length(pZ.c);
n2 = length(S.c);

% center vector
c = [pZ.c;S.c];

% generator matrix, exponent matrix and identifier vector
if isempty(pZ.G)
    if isempty(S.G)
       G = []; 
       E = []; 
       id = [];
    else
       G = [zeros(n1,size(S.G,2));S.G];
       E = S.E;
       id = S.id;
    end
else
    if isempty(S.G)
       G = [pZ.G;zeros(n2,size(pZ.G,2))];
       E = pZ.E;
       id = pZ.id;
    else
       G = blkdiag(pZ.G,S.G);
       E = blkdiag(pZ.E,S.E);
       id = [pZ.id;max(pZ.id)+S.id];
    end
end

% matrix of independent generators
GI = [];

if isempty(pZ.GI)
   if ~isempty(S.GI)
       GI = [zeros(n1,size(S.GI,2));S.GI];
   end
else
   if isempty(S.GI)
       GI = [pZ.GI;zeros(n2,size(pZ.GI,2))];
   else
       GI = blkdiag(pZ.GI,S.GI);
   end
end       
    
% generate new polyZonotope
pZ = polyZonotope(c,G,GI,E,id);

% ------------------------------ END OF CODE ------------------------------
