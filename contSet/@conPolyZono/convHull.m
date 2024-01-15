function cPZ = convHull(cPZ,S)
% convHull - Computes the convex hull of a constrained polynomial zonotope
%    and another set representation or point
%
% Syntax:
%    cPZ = convHull(cPZ)
%    cPZ = convHull(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - conPolyZono object, contSet object, or numerical vector
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
%    cPZ = conPolyZono(c,G,E,A,b,EC);
% 
%    res = convHull(cPZ);
% 
%    figure; hold on;
%    plotRandPoint(cPZ,[1,2],1000,'.k');
%    plot(cPZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linComb, enclose, polyZonotope/convHull

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    if nargin > 1
        if representsa_(cPZ,'emptySet',1e-8,'linearize',0,1)
            cPZ = S; return;
        end
        cPZ = aux_convHullMult(cPZ,S);
    else
        cPZ = aux_convHullSingle(cPZ);
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_convHullMult(cPZ,S)
% compute the convex hull of a constrained polynomial zonotope and another
% set representation

    % determine conPolyZono object
    [cPZ,S] = findClassArg(cPZ,S,'conPolyZono');
    
    % convert other set representations to constrained polynomial zonotope
    if ~isa(S,'conPolyZono')
        if isa(S,'zonotope') || isa(S,'interval') || ...
           isa(S,'polytope') || isa(S,'zonoBundle') || ...
           isa(S,'conZonotope') || isa(S,'polyZonotope') || ...
           isa(S,'capsule') || isa(S,'ellipsoid') || ...
           isa(S,'taylm')

            S = conPolyZono(S);
            
        elseif isnumeric(S)
            S = conPolyZono(S);
        else        
            % throw error for given arguments
            throw(CORAerror('CORA:noops',cPZ,S));
        end
    end

    % remove independent generators
    GI1 = cPZ.GI; GI2 = S.GI;
    cPZ.GI = zeros(dim(cPZ),0); S.GI = zeros(dim(S),0);
    
    % compute convex hull of depenent part using the linear combination
    res = aux_convHullSingle(linComb(cPZ,S));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    temp = zeros(length(cPZ.c),1);
    Z1 = zonotope(temp, GI1);
    Z2 = zonotope(temp, GI2);

    Z = enclose(Z1,Z2);

    % construct the resulting set
    res.GI = Z.G;
end


% Auxiliary Functions -----------------------------------------------------

function res = aux_convHullSingle(cPZ)
% compute the convex hull of a single constrained polynomial zonotope

    % properties
    n = dim(cPZ); a = n + 1; h = size(cPZ.G,2); p = size(cPZ.E,1);
    m = size(cPZ.A,1); q = size(cPZ.A,2);
    
    % construct auxiliary matrices
    c_ = repmat(cPZ.c,[1,a]);
    G_ = repmat(cPZ.G,[1,a]);
    temp = repmat({cPZ.E},[1,a]); E_ = blkdiag(temp{:});
    temp = repmat({ones(1,h)},[1,a]); Eh = blkdiag(temp{:});
    temp = repmat({cPZ.A},[1,a]); A_ = blkdiag(temp{:});
    b_ = repmat(cPZ.b,[a,1]);
    temp = repmat({cPZ.EC},[1,a]); R_ = blkdiag(temp{:});

    % construct resulting constrained polynomial zonotope
    c = a*cPZ.c; 
    G = [c_ G_ G_];
    E = [zeros(a*p,a),E_,E_; eye(a) zeros(a,h*a) Eh];
    A = [A_ zeros(a*m, a); zeros(1,a*q) ones(1,a)];
    b = [b_; -n];
    EC = [R_ zeros(p*a,a); zeros(a,q*a) eye(a)];
    
    % instantiate constrained polynomial zonotope
    res = conPolyZono(c,G,E,A,b,EC);

end

% ------------------------------ END OF CODE ------------------------------
