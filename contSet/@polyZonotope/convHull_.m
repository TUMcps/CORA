function S_out = convHull_(pZ,S,varargin)
% convHull_ - Computes the convex hull of a polynomial zonotope and another
%    set representation
%
% Syntax:
%    pZ = convHull_(pZ)
%    pZ = convHull_(pZ,S)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%
% Outputs:
%    S_out - convex hull
%
% Example: 
%    pZ = polyZonotope([0;0],[1 0;0 1],[],[1 3]);
%    S_out = convHull(pZ);
%
%    figure; hold on;
%    plot(S_out,[1,2],'FaceColor',[0.6 0.6 0.6],'Splits',20);
%    plot(pZ,[1,2],'FaceColor','r','Splits',6,'LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull, zonotope/enclose

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 1
    S_out = linComb(pZ,pZ);
    return
end

% ensure that numeric is second input argument
[pZ,S] = reorderNumeric(pZ,S);

% check dimensions
equalDimCheck(pZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pZ.precedence
    S_out = convHull(S,pZ,varargin{:});
    return
end

% convex hull with empty set
if representsa_(S,'emptySet',eps)
    S_out = pZ;
    return;
end

% convert to polynomial zonotope
try
    S = polyZonotope(S);
catch ME
    throw(CORAerror('CORA:noops',pZ,S));
end

% polyZonotope-polyZonotope case
S_out = aux_convHullMult(pZ,S);

end


% Auxiliary functions -----------------------------------------------------

function pZ = aux_convHullMult(pZ1,S)
% compute the convex hull of two polynomial zonotopes

    % remove independent generatros
    pZ1_ = polyZonotope(pZ1.c,pZ1.G,[],pZ1.E,pZ1.id);
    pZ2_ = polyZonotope(S.c,S.G,[],S.E,S.id);
    
    % compute convex hull of depenent part using the linear combination
    pZ = linComb(linComb(pZ1_,pZ1_),linComb(pZ2_,pZ2_));
    
    % compute convex hull of the independent part using the convex hull for
    % zonotopes
    c0 = zeros(length(pZ1.c),1);
    Z1 = zonotope(c0, pZ1.GI);
    Z2 = zonotope(c0, S.GI);

    Z = enclose(Z1,Z2);
    GI = Z.G;

    % construct the resulting set
    pZ.GI = GI;

end

% ------------------------------ END OF CODE ------------------------------
