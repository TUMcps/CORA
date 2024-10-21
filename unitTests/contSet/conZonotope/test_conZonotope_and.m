function res = test_conZonotope_and
% test_conZonotope_and - unit test function for intersection of a
%    constrained zonotope with other sets
%
% Syntax:
%    res = test_conZonotope_and
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 1: conZonotope (analytical) ----------------------------------------

% constrained zonotope 1
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono1 = conZonotope(Z,A,b);

% constrained zonotope 2
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono2 = conZonotope(Z,A,b);

% calculate intersection
intZono = cZono1 & cZono2;
V = vertices(intZono);
V = removeCollinearVertices2D(V,1e-5);

% define ground truth
V_ = [1 3 1 3;11/8 -1/8 -11/12 -7/12];

% % plot the result
% figure; hold on; axis([-3,4,-2.5,3.5]);
% plot(cZono1,[1,2],'r');
% plot(cZono2,[1,2],'b');
% plot(intZono,[1,2],'FaceColor','g');
% plot(V_(1,:),V_(2,:),'.k','MarkerSize',12);

% check correctness
assert(compareMatrices(V,V_,1e-6));

% TEST 2: polytope representing a halfspace (analytical) ------------------

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono = conZonotope(Z,A,b);

% halfspace
A = [1 -2];
b = 1;
P = polytope(A,b);

% calculate intersection
intZono = cZono & P;
V = vertices(intZono);

% define ground truth
V_ = [1 3 3 1; 0 1 2 3];

% % plot the result
% x = -4:0.1:4;
% y = (b-A(1)*x)./A(2);
% hold on
% plot(x,y,'g');
% plot(cZono,[1,2],'r');
% plot(intZono,[1,2],'b');

% check correctness
assert(compareMatrices(V,V_,1e-6));

% TEST 3: polytope representing a constrained hyperplane (analytical) -----

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono = conZonotope(Z,A,b);

% constrained hyperplane
C = [1 -2];
d = 1;
Ch = [-2 -0.5;1 0];
dh = [-4.25;2.5];
hyp = polytope(Ch,dh,C,d);

% calculate intersection
intZono = cZono & hyp;
V = vertices(intZono);
V = removeDuplicates(V,1e-6);

% define ground truth
V_ = [2 2.5;0.5 0.75];

% % plot the result
% x = -4:0.1:4;
% y = (d-C(1)*x)./C(2);
% poly = polytope([Ch;0 1],[dh;4]);
% plot(poly,[1,2],'FaceColor','m','FaceAlpha',0.5);
% hold on
% plot(x,y,'g');
% plot(cZono,[1,2],'r');
% plot(intZono,[1,2],'b','LineWidth',2);

% check correctness
assert(compareMatrices(V,V_,1e-6));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
