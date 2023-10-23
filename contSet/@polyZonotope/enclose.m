function pZ = enclose(pZ,varargin)
% enclose - encloses a polynomial zonotope and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in pZ, x2 \in pZ2, a \in [0,1] }
%    where pZ2 = M*pZ + pZplus
%
% Syntax:
%    pZ = enclose(pZ,pZ2)
%    pZ = enclose(pZ,M,pZplus)
%
% Inputs:
%    pZ - polyZonotope object
%    pZ2 - polyZonotope object
%    M - matrix for the linear transformation
%    pZplus - polyZonotope object added to the linear transformation
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    pZ2 = [1 2;-1 0]*pZ1 + [4;5];
%   
%    pZ = enclose(pZ1,pZ2);
%
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor',[0.6 0.6 0.6],'Splits',10);
%    plot(pZ1,[1,2],'FaceColor','r');
%    plot(pZ2,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    pZ2 = varargin{1};
elseif nargin == 3
    pZ2 = (varargin{1}*pZ) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% check if exponent matrices are identical
if all(size(pZ.id) == size(pZ2.id)) &&  all(pZ.id == pZ2.id) &&  ...
   all(size(pZ.E) == size(pZ2.E)) && all(all(pZ.E == pZ2.E))
    
    % compute convex hull of the dependent generators according to the
    % equation ch = (0.5 + 0.5 a)*pZ1 + (0.5 - 0.5 a)*pZ2, a \in [-1,1]
    G = [0.5 * pZ.G + 0.5 * pZ2.G, ...
         0.5 * pZ.G - 0.5 * pZ2.G, ...
         0.5 * pZ.c - 0.5 * pZ2.c]; 

    c = 0.5 * pZ.c + 0.5 * pZ2.c;

    temp = ones(1,size(pZ.E,2));
    E = [pZ.E, pZ.E; 0*temp, temp];
    E = [E, [zeros(size(E,1)-1,1); 1]];

    if ~isempty(pZ.id)
        id = [pZ.id; max(pZ.id)+1];
    else
        id = 1; 
    end

    % compute convex hull of the independent generators by using the
    % enclose function for linear zonotopes
    temp = zeros(length(pZ.c),1);
    Z1 = zonotope([temp, pZ.GI]);
    Z2 = zonotope([temp, pZ2.GI]);

    Z = enclose(Z1,Z2);
    GI = generators(Z);

    % construct resulting polynomial zonotope object
    pZ = polyZonotope(c,G,GI,E);
    pZ.id = id;
    
else
    
    % bring the exponent matrices to a common representation
    [id,E1,E2] = mergeExpMatrix(pZ.id,pZ2.id,pZ.E,pZ2.E);
    
    % extend generator end exponent matrix by center vector
    G1 = [pZ.c, pZ.G];
    G2 = [pZ2.c, pZ2.G];
    
    E1 = [zeros(length(id),1),E1];
    E2 = [zeros(length(id),1),E2];
                                      
    % compute convex hull of the dependent generators according to the
    % equation ch = (0.5 + 0.5 a)*pZ1 + (0.5 - 0.5 a)*pZ2, a \in [-1,1]
    G = 0.5 * [G1, G1, G2, -G2];
    
    h1 = size(E1,2);
    h2 = size(E2,2);
    E = [[E1, E1, E2 E2];
              zeros(1,h1), ones(1,h1), zeros(1,h2), ones(1,h2)];
          
    id = [id; max(id)+1];
    
    % compute convex hull of the independent generators by using the
    % enclose function for linear zonotopes
    temp = zeros(length(pZ.c),1);
    Z1 = zonotope([temp, pZ.GI]);
    Z2 = zonotope([temp, pZ2.GI]);

    Z = enclose(Z1,Z2);
    GI = generators(Z);
    
    % add up all generators that belong to identical exponents
    [Enew,Gnew] = removeRedundantExponents(E,G);
    
    % extract the center vector
    ind = find(sum(Enew,1) == 0);
    
    c = sum(Gnew(:,ind),2);
    Gnew(:,ind) = [];
    Enew(:,ind) = [];
    
    % construct resulting polynomial zonotope object
    pZ = polyZonotope(c,Gnew,GI,Enew);
    pZ.id = id;

end

% ------------------------------ END OF CODE ------------------------------
