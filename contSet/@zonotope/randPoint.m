function p = randPoint(Z)
% randPoint - generates a random point within a zonotope
%
% Syntax:  
%    p = randPoint(obj)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z = zonotope([1;0],rand(2,5));
%    p = randPoint(Z);
% 
%    plot(Z); hold on;
%    scatter(p(1,:),p(2,:),16,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       23-September-2008 
% Last update:   09-June-2020 (MW, remove loop)
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain number of generators
G = generators(Z);

%add generators randomly
factors = -1 + 2*rand(1,size(G,2));
p = center(Z) + sum(factors .* G, 2);


%------------- END OF CODE --------------