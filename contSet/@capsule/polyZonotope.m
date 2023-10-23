function pZ = polyZonotope(C)
% polyZonotope - enclose a capsule by a polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    C = capsule([1;2],[4;2],1.2);
%    pZ = polyZonotope(C);
% 
%    figure; hold on; box on;
%    plot(C,[1,2],'FaceColor','r');
%    plot(pZ,[1,2],'b','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, conPolyZono

% Authors:       Niklas Kochdumper
% Written:       03-October-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(C);
order = 10;

% poly. zonotope enclosure of unit ball (using spherical coordinates)
B = taylm(interval(ones(n,1)),order);

for i = 1:n-1
    if i == 1
        phi = taylm(interval(0,pi),order,['phi',num2str(i)]);
    else
        phi = taylm(interval(0,2*pi),order,['phi',num2str(i)]);
    end
    B(i) = B(i) * cos(phi);
    for j = i+1:n
       B(j) = B(j) * sin(phi); 
    end
end

B = polyZonotope(B);

T = eye(n,1);
T(1:2,1:2) = [0 1;1 0];
B = T*B;

% capsule -> convex hull of the spheres at both ends
T = eye(n); T(1,1) = -1;
d = zeros(n,1); d(1) = norm(C.g);
B1 = d + C.r*B;
B2 = -d + C.r*T*B;

pZ = linComb(B1,B2);

% consider center + direction of the generator
T = gramSchmidt(C.g);
pZ = C.c + T * pZ;
    
% ------------------------------ END OF CODE ------------------------------
