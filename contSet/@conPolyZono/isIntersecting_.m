function res = isIntersecting_(cPZ,S,type,tol,varargin)
% isIntersecting_ - determines if a constrained polynomial zonotope
%    intersects a set
%
% Syntax:
%    res = isIntersecting_(cPZ,S,type,tol)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    c = [0;0];
%    G = [8 4; 0 8];
%    E = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    EC = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    pZ1 = polyZonotope([0;0],0.5*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%
%    res1 = isIntersecting(cPZ,pZ1,'approx')
%    res2 = isIntersecting(cPZ,pZ2,'approx')
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','b','Splits',12);
%    plot(pZ1,[1,2],'g');
%    plot(pZ2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, polyZonotope/isIntersecting_, isempty

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[cPZ,S] = reorderNumeric(cPZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cPZ.precedence
    res = isIntersecting_(S,cPZ,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(cPZ,S,type,tol,0,false,false);
    return
end

% currently, no exact algorithm
if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',cPZ,S));
end

% call function for other set representations
if isa(S,'conPolyZono') || isa(S,'polyZonotope') || isa(S,'capsule')
    res = aux_isIntersecting_cPZ_pZ_C(cPZ,S,type,tol);
    return
end

% convert to polytope
if isa(S,'polytope')
    res = aux_isIntersecting_polytope(cPZ,S,type,tol);
    return
end
if isa(S,'interval')
    res = aux_isIntersecting_polytope(cPZ,polytope(S),type,tol);
    return
end

% enclose by constrained zonotope
if isa(S,'zonotope') || isa(S,'zonoBundle') || isa(S,'conZonotope')
    res = isIntersecting_(conZonotope(cPZ),conZonotope(S),type,tol);
    return
end

throw(CORAerror('CORA:noops',cPZ,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isIntersecting_cPZ_pZ_C(cPZ,S,type,tol)

% fast check based on zonotope enclosure
if ~isIntersecting_(zonotope(cPZ),zonotope(S),type,tol)
    res = false;
    return
end

% convert second set to constrained polynomial zonotope
S = conPolyZono(S);

% compute intersection of the two sets
I = and_(cPZ,S,'exact');

% check if the intersection is empty
res = ~representsa_(I,'emptySet',tol);

end

function res = aux_isIntersecting_polytope(cPZ,P,type,tol)

res = true;

% check intersection with hyperplane
if representsa_(P,'conHyperplane',tol)
    I = supportFunc_(cPZ,P.Ae_.val(1,:)','range','interval',8,1e-3);
    if ~contains_(I,P.be_.val(1),'exact',tol,0,false,false)
        res = false;
        return
    end
end

% loop over all halfspaces
for i = 1:size(P.A,1)
    if P.b(i) < supportFunc_(cPZ,P.A(i,:)','lower','interval',8,1e-3)
        res = false;
        return
    end
end

end

% ------------------------------ END OF CODE ------------------------------
