function pZ = polyZonotope(E)
% polyZonotope - enclose an ellipsoid by a polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    E = ellipsoid([4 2;2 4],[1;1]);
%    pZ = polyZonotope(E);
% 
%    figure; hold on; xlim([-2,4]); ylim([-2,4]);
%    plot(E,[1,2],'FaceColor','r');
%    plot(pZ,[1,2],'b');
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

% read out dimension, fix order
n = dim(E);
order = 10;

% poly. zonotope enclosure of unit ball (using spherical coordinates)
r = taylm(interval(0,1),order,'r');
B = r * ones(n,1);

for i=1:n-1
    phi = taylm(interval(0,2*pi),order,['phi',num2str(i)]);
    B(i) = B(i) * cos(phi);
    for j = i+1:n
        B(j) = B(j) * sin(phi); 
    end
end

% instantiate polynomial zonotope
B = polyZonotope(B);

% ellipsoid -> linear transformation of the unit ball
[V,D] = eig(E.Q);
pZ = E.q + (sqrt(D)*V)'*B;
    
% ------------------------------ END OF CODE ------------------------------
