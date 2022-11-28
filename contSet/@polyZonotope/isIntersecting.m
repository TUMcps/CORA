function res = isIntersecting(pZ,S,varargin)
% isIntersecting - determines if a polynomial zonotope intersects a set
%
% Syntax:  
%    res = isIntersecting(pZ,S)
%    res = isIntersecting(pZ,S,type)
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
% See also: conZonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_isIntersecting('polyZonotope',pZ,S,varargin{:});

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign values
    pZ = vars{1}; S = vars{2}; type = vars{3};
end


if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',S,pZ));
end

% get polyZonotope object
if ~isa(pZ,'polyZonotope')
    temp = pZ;
    pZ = S;
    S = temp;
end

% call function for other set representations
if isa(S,'polyZonotope')
    
    % fast check based on zonotope enclosure
    res = isIntersecting(zonotope(pZ),zonotope(S));
    
    if ~res
        return; 
    end
    
    % construct polynomial constraint for the intersection
    c = pZ.c - S.c;
    G = [pZ.G -S.G];
    Grest = [pZ.Grest,-S.Grest];
    expMat = blkdiag(pZ.expMat,S.expMat);
    
    % contract the factor domain \alpha_k in [-1,1] based on the
    % polynomial constraint
    n = size(expMat,1) + size(Grest,2);
    dom = interval(-ones(n,1),ones(n,1));
    
    int = contractPoly(c,G,Grest,expMat,dom,'linearize',3,7);
    
    % check if polynomial zonotopes are intersecting
    if isempty(int)
       res = false;
    else
       res = true;
    end
    
elseif isa(S,'halfspace') || isa(S,'conHyperplane') || ...
       isa(S,'mptPolytope') || isa(S,'zonotope') || ...
       isa(S,'interval') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope') || isa(S,'ellipsoid') || ...
       isa(S,'conPolyZono')

    res = isIntersecting(S,pZ,type);

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',pZ,S));
end

%------------- END OF CODE --------------
