function p = randPointExtreme(Z)
% randPointExtreme - generates a random extreme point of a zonotope
%
% Syntax:  
%    p = randPointExtreme(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    zono = zonotope.generateRandom(2,[],3);
%    p = randPointExtreme(zono);
%
%    figure
%    hold on
%    plot(zono,[1,2],'r');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      14-May-2009
% Last update:  09-June-2020 (MW, remove loop)
% Last revision:---

%------------- BEGIN CODE --------------

% obtain generator matrix
G = generators(Z);

% add generators with random factors -1 or 1
factors = sign(-1 + 2*rand(1,size(G,2)));
p = center(Z) + sum(factors .* G, 2);

%------------- END OF CODE --------------