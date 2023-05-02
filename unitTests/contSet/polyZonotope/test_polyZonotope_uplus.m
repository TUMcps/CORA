function res = test_polyZonotope_uplus
% test_polyZonotope_uplus - unit test function of uminus
%
% Syntax:  
%    res = test_polyZonotope_uplus
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Tobias Ladner
% Written:      06-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = true(0);

% init
c = [1;2];
G = [1 0 1; 0 1 1];
Grest = [1;2];
E = [1 0 4; 2 1 3];
pZ = polyZonotope(c,G,Grest,E);

% plus
ppZ = +pZ;
resvec(end+1) = all(ppZ.c == c, 'all');
resvec(end+1) = all(ppZ.G == G, 'all');
resvec(end+1) = all(ppZ.Grest == Grest, 'all');

% compare with pZ
resvec(end+1) = isequal(ppZ, pZ);

% test empty case
resvec(end+1) = isemptyobject(+polyZonotope());

% add results
res = all(resvec);

%------------- END OF CODE --------------