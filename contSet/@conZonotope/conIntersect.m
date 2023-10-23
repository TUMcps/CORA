function res = conIntersect(cZ1,cZ2,M)
% conIntersect - Adds the constraint that the linear transformation M*cZ1
%    of a constrained zonotope cZ1 has to intersect the constrained 
%    zonotope cZ2 (see Eq. (13) in [1] with cZ1=Z, cZ2=Y, M=R)
%
% Syntax:
%    res = conIntersect(cZ1,cZ2,M)
%
% Inputs:
%    cZ1 - conZonotope object
%    cZ2 - conZonotope object
%    M - transformation matrix for the linear transformation M*cZ1
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    cZ1 = conZonotope([0 1.5 -1.5 0.5;0 1 0.5 -1],[1 1 1],-1);
%    cZ2 = conZonotope([6 2 0;-2 0 1]);
%    M = [2 1;-1 0];
% 
%    cZ_ = conIntersect(cZ1,cZ2,M);
% 
%    figure; hold on;
%    plot(cZ1); plot(cZ2);
%    plot(cZ_,[1,2],'r--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and, linearSysDT/observe
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
       
% get object properties
cz = cZ1.c; cy = cZ2.c; Gz = cZ1.G; Gy = cZ2.G;

% construct resulting set according to Eq. (13) in [1]
c = cz;
G = [Gz, zeros(size(Gz,1),size(Gy,2))];

if ~isempty(cZ1.A)
   if ~isempty(cZ2.A)
       A = [blkdiag(cZ1.A,cZ2.A);[M*Gz -Gy]];
       b = [cZ1.b; cZ2.b; cy - M*cz];
   else
       A = [[cZ1.A,zeros(size(cZ1.A,1),size(Gy,2))];[M*Gz -Gy]];
       b = [cZ1.b; cy - M*cz];
   end
else
   if ~isempty(cZ2.A)
       A = [[zeros(size(cZ2.A,1),size(Gz,2)) cZ2.A];[M*Gz -Gy]];
       b = [cZ2.b; cy - M*cz];
   else
       A = [M*Gz -Gy];
       b = cy - M*cz;
   end
end

% instantiate intersection
res = conZonotope(c,G,A,b);

% ------------------------------ END OF CODE ------------------------------
