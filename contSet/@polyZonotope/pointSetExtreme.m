function points = pointSetExtreme(pZ)
% pointSet - Computes the set of all extreme points (all factors 1 or -1)
%            of a polynomial zonotope
%
% Syntax:  
%    points = pointSetExtreme(pZ)
%
% Inputs:
%    pZ - polyZonotope object (n dimensional)
%
% Outputs:
%    points - extreme point set (dimension: [n,M]) 
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1 1;0 2 1 2],[],[1 0 3 1;0 1 0 2]);
%    points = pointSetExtreme(pZ);
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: pointSet

% Author:       Niklas Kochdumper
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% determine the number of factors
N = size(pZ.expMat,1);
M = size(pZ.Grest,2);
P = N+M;

% determine all possible extreme points (all factors equal to 1 or -1)
temp = interval(-ones(P,1),ones(P,1));
V = vertices(temp);

points = zeros(size(pZ.G,1),size(V,2));

% loop over all extreme points
for i = 1:size(V,2)

    % Part 1: dependent generators
    if ~isempty(pZ.G)
        p = pZ.c;
        beta = V(1:N,i);
        fact = prod(beta.^pZ.expMat,1);
        p = p + pZ.G * fact';
    end

    % Part 2: independent generators
    if ~isempty(pZ.Grest)
        alpha = V(N+1:end,i);
        p = p + pZ.Grest * alpha;
    end

    points(:,i) = p;
end

%------------- END OF CODE --------------