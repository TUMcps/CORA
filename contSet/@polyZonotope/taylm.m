function res = taylm(pZ)
% taylm - encloses a polynomial zonotope with a Taylor model
%
% Syntax:
%    res = taylm(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - taylm object
%
% Example:
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    tay = taylm(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope

% Authors:       Niklas Kochdumper
% Written:       13-August-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% introduce independent factors as new dependent factors
G = [pZ.G,pZ.GI];
E = blkdiag(pZ.E,eye(size(pZ.GI,2)));

% create Taylor model for factors
p = length(pZ.id) + size(pZ.GI,2);
int = interval(-ones(p,1),ones(p,1));

tay = taylm(int);

% convert polyZonotope object to taylor model
res = pZ.c;

for i = 1:size(G,2)
    temp = 1;
    for j = 1:size(E,1)
       temp = temp * tay(j)^E(j,i); 
    end
    res = res + G(:,i) * temp;
end

% ------------------------------ END OF CODE ------------------------------
