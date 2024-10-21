function S_out = convHull_(cPZ,S,varargin)
% convHull_ - Computes the convex hull of a constrained polynomial zonotope
%    and another set representation or point
%
% Syntax:
%    cPZ = convHull_(cPZ)
%    cPZ = convHull_(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - conPolyZono object, contSet object, or numerical vector
%
% Outputs:
%    S_out - convex hull
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
% See also: contSet/convHull, linComb, enclose, polyZonotope/convHull_

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% convex hull of set itself
if nargin == 1
    S_out = aux_convHullSingle(cPZ);
    return
end

% ensure that numeric is second input argument
[cPZ,S] = reorderNumeric(cPZ,S);

% check dimensions
equalDimCheck(cPZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cPZ.precedence
    S_out = convHull(S,cPZ,varargin{:});
    return
end

% check for emptiness
if representsa_(cPZ,'emptySet',1e-8,'linearize',0,1)
    S_out = S;
    return
end

% convex hull of constrained polynomial zonotope and other set: always via
% conversion of second set to conPolyZono object
try
    S = conPolyZono(S);
catch ME
    throw(CORAerror('CORA:noops',cPZ,S));
end
S_out = aux_convHullMult(cPZ,S);

end


% Auxiliary functions -----------------------------------------------------

function S_out = aux_convHullMult(cPZ,S)
% compute the convex hull of a constrained polynomial zonotope and another
% set representation

    % read out dimension
    n = dim(cPZ);

    % remove independent generators
    GI1 = cPZ.GI; GI2 = S.GI;
    cPZ.GI = zeros(n,0); S.GI = zeros(n,0);
    
    % compute convex hull of depenent part using the linear combination
    S_out = aux_convHullSingle(linComb(cPZ,S));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    Z = enclose(zonotope(zeros(n,1),GI1), zonotope(zeros(n,1),GI2));

    % construct the resulting set
    S_out.GI = Z.G;
end


% Auxiliary Functions -----------------------------------------------------

function S_out = aux_convHullSingle(cPZ)
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
    S_out = conPolyZono(c,G,E,A,b,EC);

end

% ------------------------------ END OF CODE ------------------------------
