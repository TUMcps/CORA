function res = linComb(cPZ,S)
% linComb - Computes the linear combination of a constrained polynomial
%    zonotope and another set representation or point
%
% Syntax:
%    res = linComb(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object or numerical vector
%
% Outputs:
%    res - conPolyZono object enclosing cPZ1 and cPZ2
%
% Example: 
%    c = [0;0];
%    G = [1 1 0;0 0 1];
%    E = [1 0 0; 0 1 0; 0 0 1; 0 0 0;0 0 0];
%    A = [1 1 -0.5 0; 0 0 0 1];
%    b = [0.5; 1];
%    EC = [0 0 0 1; 2 0 0 0; 0 2 0 0; 0 0 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    Z = zonotope([0;-2],[0.5 0.5; 0.5 -0.5]);
%   
%    res = linComb(cPZ,Z);
%
%    figure; hold on;
%    plot(res,[1,2],'FaceColor','r','Splits',12);
%    plot(cPZ,[1,2],'FaceColor','b','Splits',12);
%    plot(Z,[1,2],'FaceColor','g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: convHull, enclose, polyZonotope/linComb

% Authors:       Niklas Kochdumper
% Written:       21-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine conPolyZono object
if ~isa(cPZ,'conPolyZono')
    temp = cPZ;
    cPZ = S;
    S = temp;
end

% convert other set representations to constrained polynomial zonotope
if ~isa(S,'conPolyZono')
    if isa(S,'zonotope') || isa(S,'interval') || ...
       isa(S,'polytope') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope') || isa(S,'polyZonotope') || ...
       isa(S,'capsule') || isa(S,'ellipsoid') || ...
       isa(S,'taylm')

        S = conPolyZono(S);
        
    elseif isnumeric(S)
        S = conPolyZonotope(S,[],[]);
    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',cPZ,S));
    end
end

% call linComb method for polynomial zonotopes
pZ1 = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
pZ2 = polyZonotope(S.c,S.G,S.GI,S.E,S.id);

pZ = linComb(pZ1,pZ2);

% convert to constrained polynomial zonotope
res = conPolyZono(pZ);

% update constraints
res.A = blkdiag(cPZ.A,S.A);
res.b = [cPZ.b;S.b];

if isempty(cPZ.A)
    if ~isempty(S.A)
        temp = zeros(length(cPZ.id),size(S.EC,2));
        EC = [temp;S.EC];
        res.EC = [EC; zeros(1,size(res.A,2))];
    end
else
    if isempty(S.A)
        temp = zeros(length(S.id),size(cPZ.EC,2));
        EC = [cPZ.EC;temp];
        res.EC = [EC; zeros(1,size(res.A,2))];
    else
        EC = blkdiag(cPZ.EC,S.EC);
        res.EC = [EC; zeros(1,size(res.A,2))];
    end
end

% ------------------------------ END OF CODE ------------------------------
