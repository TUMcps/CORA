function V = underapproximate(Z,S)
% underapproximate - returns the vertices of an underapproximation. The
%    underapproximation is computed by finding the vertices that are
%    extreme in the direction of a set of vectors, stored in the matrix S.
%    If S is not specified, it is constructed by the vectors spanning an
%    over-approximative parallelotope.
%
% Syntax:
%    V = underapproximate(Z,S)
%
% Inputs:
%    Z - zonotope object
%    S - matrix of direction vectors
%
% Outputs:
%    V - vertices
%
% Example: ---
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices

% Authors:       Matthias Althoff
% Written:       19-July-2010
% Last update:   28-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract center and generators
c = Z.c;
G = Z.G;

if nargin==1
    [rows,cols] = size(G);
    if rows>=cols
        S = G;
    else
        S = dominantDirections(Z);
    end
end

%obtain extreme vertices along directions in S
V=[0*S, 0*S];
for i=1:length(S(1,:))
    posVertex = c;
    negVertex = c;
    for iGen=1:length(G(1,:))
        s = sign(S(:,i).'*G(:,iGen));
        posVertex = posVertex + s*G(:,iGen);
        negVertex = negVertex - s*G(:,iGen);
    end
    V(:,2*i-1) = posVertex;
    V(:,2*i) = negVertex;
end

% ------------------------------ END OF CODE ------------------------------
