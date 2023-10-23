function Zres = intersectStripAlamo_unitTest(Z,C,phi,d)
% intersectStripAlamo_unitTest - computes the intersection between a zonotope and
% a strip defined as |C x-d| <= phi according to [1]
%
% Syntax:
%    Zres = intersectStripAlamo_unitTest(Z,C,phi,d)
%
% Inputs:
%    Z - zonotope object
%    C - normal vector of strip
%    phi - width of strip
%    d - center of strip
%
% Outputs:
%    res - true/false
%
% Example:
%    % strip and zonotope
%    C = [1 0];
%    phi = 5;
%    d = -2;
% 
%    Z = zonotope([1 2 2 2 6 2 8; 1 2 2 0 5 0 6 ]);
%    res_zono = intersectStripAlamo(Z,C,phi,d);
% 
%    % just for comparison
%    poly = polytope([1 0;-1 0],[3;7]);
%    Zpoly = Z & poly;
% 
%    figure; hold on 
%    plot(Z,[1 2],'r-+');
%    plot(poly,[1 2],'r-*');
%    plot(Zpoly,[1 2],'b-+');
%    plot(res_zono,[1 2],'b-*');
% 
%    legend('zonotope','strip','zono&poly','zonoStrips');
%
% References:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035-1043,
%        2005.

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract generators
G = generators(Z);

% auxiliary variables
aux1 = G*G'; 
aux2 = aux1*C';
aux3 = C*aux1*C' + phi^2; 

% return lambda
lambda = aux2/aux3;

%% resulting zonotope
% extract center
c = center(Z);
% new center
c_new = c + lambda*(d - C*c);

% new generators
I = eye(length(c));
part1 = (I - lambda*C)*G;
part2 = phi*lambda;
G_new = [part1 part2];

% resulting zonotope
Zres = zonotope([c_new G_new]);

% ------------------------------ END OF CODE ------------------------------
