function cZ = intervalMultiplication(cZ,I)
% intervalMultiplication - computes the multiplication of an interval with
%    a constrained zonotope
%
% Syntax:
%    res = intervalMultiplication(obj,I)
%
% Inputs:
%    obj - conZonotope object 
%    I - interval object
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    I = interval([0 1; 1 0], [1 2; 2 1]);
%    cZres = I * cZ;
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(cZres,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, interval/mtimes, zonotope/intervalMultiplication

% Authors:       Niklas Kochdumper
% Written:       03-July-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% center and radius of interval matrix
m = center(I);
r = rad(I);

% absolute value of zonotope center and generators
Zabssum = sum(abs([cZ.c,cZ.G]),2);

% construct resulting conZonotope object
cZ.c = m*cZ.c;
cZ.G = [m*cZ.G, diag(r*Zabssum)];
cZ.A = [cZ.A, zeros(size(cZ.A,1),size(Zabssum,1))];

cZ.ksi = [];
cZ.R = [];

% ------------------------------ END OF CODE ------------------------------
