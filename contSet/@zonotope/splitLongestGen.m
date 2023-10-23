function Znew = splitLongestGen(Z)
% splitLongestGen - splits the longest generator
%
% Syntax:
%    Znew = splitLongestGen(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Znew - cell array storing the split zonotope objects
%
% Example: 
%
%    Z = zonotope([1;0],[1 3 -2 -1; 0 2 -1 1]);
%    Znew = splitLongestGen(Z);
%
%    figure; hold on; box on;
%    plot(Z);
%    plot(Znew{1},[1,2],'--r');
%    plot(Znew{2},[1,2],'--g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/splitFirstGen, zonotope/split

% Authors:       Niklas Kochdumper
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% object properties
c = center(Z);
G = generators(Z);

% determine longest generator   
len = sum(G.^2, 1);
[~,ind] = max(len);

% split longest generator
c1 = c + 0.5*G(:,ind);
G1 = G;
G1(:,ind) = 0.5*G1(:,ind);

c2 = c - 0.5*G(:,ind);
G2 = G;
G2(:,ind) = 0.5*G2(:,ind);

Znew = {zonotope(c1,G1), zonotope(c2,G2)};

% ------------------------------ END OF CODE ------------------------------
