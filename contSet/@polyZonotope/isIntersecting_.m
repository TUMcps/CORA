function res = isIntersecting_(pZ,S,type,tol,varargin)
% isIntersecting_ - determines if a polynomial zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(pZ,S,type,tol)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    pZ2 = polyZonotope([1;4],[1 1 1 1;1 0 -1 1],[],[1 0 1 2;0 1 1 0]);
%    P1 = polytope([-1 1],-2);
%    P2 = polytope([-1 1],-4);
%
%    isIntersecting(pZ1,pZ2,'approx')
%    isIntersecting(pZ1,P1,'approx')
%    isIntersecting(pZ1,P2,'approx')
%
%    figure; hold on;
%    plot(pZ1,[1,2],'b');
%    plot(pZ2,[1,2],'r');
%
%    figure; hold on; xlim([-4,6]); ylim([-5,5]);
%    plot(P1,[1,2],'b');
%    plot(pZ1,[1,2],'FaceColor','g');
%
%    figure; hold on; xlim([-4,6]); ylim([-5,5]);
%    plot(P2,[1,2],'b');
%    plot(pZ1,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, conZonotope/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       21-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[pZ,S] = reorderNumeric(pZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pZ.precedence
    res = isIntersecting_(S,pZ,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(pZ,S,type,tol);
    return
end

% sets must not be empty
if representsa_(pZ,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% check if polynomial zonotope actually represents a zonotope
if representsa_(pZ,'zonotope',tol)
    res = isIntersecting_(zonotope(pZ),S,type,tol);
    return
end

% currently, no exact algorithm supported
if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',pZ,S));
end

% polyZonotope intersection check
if isa(S,'polyZonotope')
    res = aux_isIntersecting_polyZonotope(pZ,S,type,tol);
    return
end

% convert to polytope if easily possible
if isa(S,'polytope')
    res = aux_isIntersecting_polytope(pZ,S,type,tol);
    return
end
if isa(S,'interval') 
    res = aux_isIntersecting_polytope(pZ,polytope(S),tol);
    return
end

throw(CORAerror('CORA:noops',pZ,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isIntersecting_polytope(pZ,P,type,tol)

res = true;

% check intersection with hyperplane
if representsa_(P,'conHyperplane',tol)
    I = supportFunc_(pZ,P.Ae(1,:)','range','interval',8,1e-3);
    if ~contains_(I,P.be(1),'exact',tol)
        res = false;
        return
    end
end

% loop over all halfspaces
for i = 1:size(P.A,1)
    if P.b(i) < supportFunc_(pZ,P.A(i,:)','lower','interval',8,1e-3)
        res = false;
        return
    end
end

end

function res = aux_isIntersecting_polyZonotope(pZ,S,type,tol)

% fast check based on zonotope enclosure
if ~isIntersecting_(zonotope(pZ),zonotope(S),type,tol)
    res = false;
    return
end

% construct polynomial constraint for the intersection
c = pZ.c - S.c;
G = [pZ.G -S.G];
GI = [pZ.GI,-S.GI];
E = blkdiag(pZ.E,S.E);

% contract the factor domain \alpha_k in [-1,1] based on the
% polynomial constraint
n = size(E,1) + size(GI,2);
dom = interval(-ones(n,1),ones(n,1));

I = contractPoly(c,G,GI,E,dom,'linearize',3,7);

% check if polynomial zonotopes are intersecting
res = ~representsa_(I,'emptySet',tol);

end

% ------------------------------ END OF CODE ------------------------------
