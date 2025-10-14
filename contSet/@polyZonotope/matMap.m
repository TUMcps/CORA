function pZ = matMap(pZ1,pZ2,n,k,m)
% matMap - computes the matrix multiplication of the two given matrix 
%    polynomial zonotopes [1, Lemma 1]
%
% Syntax:
%    pZ = matMap(pZ1,pZ2,n,k,m)
%
% Inputs:
%    pZ1 - polyZonotope object (representing a n x k matrix)
%    pZ2 - polyZonotope object (representing a k x m matrix)
%    n,k,m - dimensions
%
% Outputs:
%    pZ - polyZonotope object (representing a n x m matrix)
%
% Example:
%    n = 2; k = 3; m = 2;
%    M1 = randn(n,k); M2 = randn(k,m);
%    epsilon = 0.01;
%    pZ1 = polyZonotope(reshape(M1,n*k,1), eye(n*k)*epsilon);
%    pZ2 = polyZonotope(reshape(M2,k*m,1), eye(k*m)*epsilon);
%    pZ_res = matMap(pZ1, pZ2, n, k, m)
%
% Reference:
%    [1] Ladner et al. "Formal Verification of Graph Convolutional Networks
%        with Uncertain Node Features and Uncertain Graph Structure".
%        TMLR. 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap

% Authors:       Tobias Ladner
% Written:       07-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(5,5);
aux_validateInput(pZ1,pZ2,n,k,m)

% read out properties
c1 = pZ1.c; G1 = pZ1.G; GI1 = pZ1.GI; E1 = pZ1.E; id1 = pZ1.id;
c2 = pZ2.c; G2 = pZ2.G; GI2 = pZ2.GI; E2 = pZ2.E; id2 = pZ2.id;

% prepare for multiplication
h1 = size(G1,2); h2 = size(G2,2); hi1 = size(GI1,2); hi2 = size(GI2,2);
c1 = reshape(c1,n,k); G1 = reshape(G1,n,k,h1,1); GI1 = reshape(GI1,n,k,hi1,1);
c2 = reshape(c2,k,m); G2 = reshape(G2,k,m,1,h2); GI2 = reshape(GI2,k,m,1,hi2);

% compute output center ---
c1c2 = c1 * c2;
c1c2 = reshape(c1c2,n*m,1);
c = c1c2;

% compute output generator matrix ---
G1c2 = pagemtimes(G1, c2);
c1G2 = pagemtimes(c1, G2);
G1G2 = pagemtimes(G1, G2);

G1c2 = reshape(G1c2,n*m,h1);
c1G2 = reshape(c1G2,n*m,h2);
G1G2 = reshape(G1G2,n*m,h1*h2);

G = [G1c2, c1G2, G1G2];

% compute output exponent matrix ---
[id,E1,E2] = mergeExpMatrix(id1,id2,E1,E2);
p = numel(id);
E = [E1, E2];
if ~isempty(E1) && ~isempty(E2)
    % sum each column of E1 to E2 and concatenate
    E = [E, reshape(E1 + reshape(E2,p,1,h2),p,h1*h2)];
end

% compute output independent generator matrix ---
GI1c2 = pagemtimes(GI1, c2);
c1GI2 = pagemtimes(c1, GI2);
GI1GI2 = pagemtimes(GI1, GI2);
G1GI2 = pagemtimes(G1,GI2);
GI1G2 = pagemtimes(GI1,G2);

GI1c2 = reshape(GI1c2,n*m,hi1);
c1GI2 = reshape(c1GI2,n*m,hi2);
GI1GI2 = reshape(GI1GI2,n*m,hi1*hi2);
G1GI2 = reshape(G1GI2,n*m,h1*hi2);
GI1G2 = reshape(GI1G2,n*m,hi1*h2);

GI = [GI1c2, c1GI2, GI1GI2, G1GI2, GI1G2];

% init output polynomial zonotope ---
pZ = polyZonotope(c, G, GI, E, id);
 
end


% Auxiliary functions -----------------------------------------------------

function aux_validateInput(pZ1,pZ2,n,k,m)
    % validate inputs
    inputArgsCheck({ ...
        {pZ1,'att','polyZonotope'}, ...
        {pZ2,'att','polyZonotope'}, ...
        {n,'att','numeric',{'scalar','integer'}}, ...
        {k,'att','numeric',{'scalar','integer'}}, ...
        {m,'att','numeric',{'scalar','integer'}}, ...
    })
    % check dimensions
    if dim(pZ1) ~= n*k
        throw(CORAerror('CORA:wrongValue','third/fourth','Dimension of pZ1 must match given n*k.'))
    end
    if dim(pZ2) ~= k*m
        throw(CORAerror('CORA:wrongValue','fourth/fifth','Dimension of pZ2 must match given k*m.'))
    end
end

% ------------------------------ END OF CODE ------------------------------
