function res = isIntersecting_(pZ,S,type,varargin)
% isIntersecting_ - determines if a polynomial zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(pZ,S)
%    res = isIntersecting_(pZ,S,type)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    pZ2 = polyZonotope([1;4],[1 1 1 1;1 0 -1 1],[],[1 0 1 2;0 1 1 0]);
%    hs1 = halfspace([-1 1],-2);
%    hs2 = halfspace([-1 1],-4);
%
%    isIntersecting(pZ1,pZ2,'approx')
%    isIntersecting(pZ1,hs1,'approx')
%    isIntersecting(pZ1,hs2,'approx')
%
%    figure; hold on;
%    plot(pZ1,[1,2],'b');
%    plot(pZ2,[1,2],'r');
%
%    figure; hold on; xlim([-4,6]); ylim([-5,5]);
%    plot(hs1,[1,2],'b');
%    plot(pZ1,[1,2],'FaceColor','g');
%
%    figure; hold on; xlim([-4,6]); ylim([-5,5]);
%    plot(hs2,[1,2],'b');
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

if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',S,pZ));
end

% get polyZonotope object
[pZ,S] = findClassArg(pZ,S,'polyZonotope');

% call function for other set representations
if isa(S,'polyZonotope')
    
    % fast check based on zonotope enclosure
    if ~isIntersecting_(zonotope(pZ),zonotope(S),type)
        res = false; return
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
    res = ~representsa_(I,'emptySet',eps);
    
elseif isa(S,'halfspace') || isa(S,'conHyperplane') || ...
       isa(S,'polytope') || isa(S,'zonotope') || ...
       isa(S,'interval') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope') || isa(S,'ellipsoid') || ...
       isa(S,'conPolyZono')

    res = isIntersecting_(S,pZ,type);

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',pZ,S));
end

% ------------------------------ END OF CODE ------------------------------
