function pZ = linComb(pZ1,S)
% linComb - computes the linear combination of two polynomial zonotopes
%
% Syntax:
%    pZ = linComb(pZ1,S)
%
% Inputs:
%    pZ1 - polyZonotope object
%    S - contSet object
%
% Outputs:
%    pZ - polyZonotope enclosing pZ1 and pZ2
%
% Example: 
%    pZ1 = polyZonotope([-2;-2],[2 0 1;0 2 1],[],[1 0 3;0 1 1]);
%    pZ2 = polyZonotope([3;3],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%   
%    pZ = linComb(pZ1,pZ2);
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
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine polyZonotope object
if ~isa(pZ1,'polyZonotope')
    temp = pZ1;
    pZ1 = S;
    S = temp;
end

% convert other set representations to polynomial zonotopes
if ~isa(S,'polyZonotope')
    if isa(S,'zonotope') || isa(S,'interval') || ...
       isa(S,'polytope') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope')

        S = polyZonotope(S);
        
    elseif isa(S,'conPolyZono')
        
        pZ = linComb(S,pZ1);
        return;
        
    elseif isnumeric(S)
        
        S = polyZonotope(S,[],[],[]);
        
    else        
        % throw error for given arguments
        throw(CORAerror('CORA:noops',pZ1,S));
    end
end

% extend generator end exponent matrix by center vector
G1 = [pZ1.c, pZ1.G];
G2 = [S.c, S.G];

E1 = [zeros(size(pZ1.E,1),1),pZ1.E];
E2 = [zeros(size(S.E,1),1),S.E];

% compute linear comb. of the dependent generators according to the
% equation comb = (0.5 + 0.5 a)*pZ1 + (0.5 - 0.5 a)*pZ2, a \in [-1,1]
G = 0.5 * [G1, G1, G2, -G2];

h1 = size(E1,2);
h2 = size(E2,2);

zero1 = zeros(size(E1,1),size(E2,2));
zero2 = zeros(size(E2,1),size(E1,2));

E = [E1, E1, zero1, zero1; ...
          zero2, zero2, E2 E2; ...
          zeros(1,h1), ones(1,h1), zeros(1,h2), ones(1,h2)];

% compute convex hull of the independent generators by using the
% enclose function for linear zonotopes
temp = zeros(length(pZ1.c),1);
Z1 = zonotope([temp, pZ1.GI]);
Z2 = zonotope([temp, S.GI]);

Z = enclose(Z1,Z2);
GI = Z.G;

% add up all generators that belong to identical exponents
[Enew,Gnew] = removeRedundantExponents(E,G);

% extract the center vector
ind = find(sum(Enew,1) == 0);

c = sum(Gnew(:,ind),2);
Gnew(:,ind) = [];
Enew(:,ind) = [];

% construct resulting polynomial zonotope object
pZ = polyZonotope(c,Gnew,GI,Enew);
pZ.id = (1:size(Enew,1))';

% ------------------------------ END OF CODE ------------------------------
