function [pZ] = randProbZonotope(dim,detGenerators,probGenerators)
% randProbZonotope - Generates a random probabilistic zonotope.
%
% Syntax:  
%    [Obj] = randProbZonotope(dim,detGenerators,probGenerators)
%
% Inputs:
%    dim - dimension of the vector space
%    detGenerators - number of deterministic generators
%    probGenerators - number of probabilistic generators
%
% Outputs:
%    Obj - random probabilistic zonotope object structure
%
% Example: 
%    pZ = probZonotope(2,4,3);
%    plot(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 23-August-2007
% Last update: 04-March-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%generators
Z=rand(dim,detGenerators+1);
g=-1+2*rand(dim,probGenerators);
gamma=3;

pZ=probZonotope(Z,g,gamma);

%------------- END OF CODE --------------