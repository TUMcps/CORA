function completed = example_setRepresentationsAndOperations()
% updated: 20-March-2015, MA
% updated: 30-August-2016, MA

%zonotope------------------------------------------------------------------
Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2
A = [0.5 1; 1 0.5]; % numerical matrix A

Z3 = Z1 + Z2; % Minkowski addition
Z4 = A*Z3; % linear map

figure; hold on
plot(Z1,[1 2],'b'); % plot Z1 in blue
plot(Z2,[1 2],'g'); % plot Z2 in green
plot(Z3,[1 2],'r'); % plot Z3 in red
plot(Z4,[1 2],'k'); % plot Z4 in black

P = polytope(Z4) % convert to and display halfspace representation
IH = interval(Z4) % convert to and display interval

figure; hold on
plot(Z4); % plot Z4
plot(IH,[1 2],'g'); % plot intervalhull in green
%--------------------------------------------------------------------------

%zonotopeBundle------------------------------------------------------------
Z{1} = Z1;
Z{2} = Z2;
Zb = zonotopeBundle(Z); % instantiate zonotope bundle from Z1, Z2
vol = volume(Zb) % compute and display volume of zonotope bundle

figure; hold on
plot(Z1); % plot Z1 
plot(Z2); % plot Z2 
plotFilled(Zb,[1 2],[.675 .675 .675],'EdgeColor','none'); % plot Zb in gray
%--------------------------------------------------------------------------

%quadZonotope--------------------------------------------------------------
c = [0;0]; % starting point
E1 = diag([-1,0.5]); % generators of factors with identical indices
E2 = [1 1; 0.5 0.3]; % generators of factors with identical indices
F = [-0.5; 1]; % generators of factors with different indices
G = [0.3; 0.3]; % independent generators

qZ = quadZonotope(c,E1,E2,F,G); % instantiate quadratic zonotope
Z = zonotope(qZ) % over-approximate by a zonotope

figure; hold on
plot(Z); % plot Z 
plotFilled(qZ,[1 2],7,[],[.6 .6 .6],'EdgeColor','none'); % plot qZ 
%--------------------------------------------------------------------------

%probZonotope--------------------------------------------------------------
Z1=[10 ; 0]; % uncertain center
Z2=[0.6 1.2  ; 0.6 -1.2]; % generators with normally distributed factors
pZ=probZonotope(Z1,Z2,2); % probabilistic zonotope

M=[-1 -1;1 -1]*0.2; % mapping matrix
pZencl = enclose(pZ,M); % probabilistic enclosur of pZ and M*pZ

figure('renderer','zbuffer')
hold on
plot(pZ,'dark'); % plot pZ
plot(expm(M)*pZ,'light'); % plot expm(M)*pZ
plot(pZencl,'mesh') % plot enclosure

campos([-3,-51,1]); %set camera position
drawnow; % draw 3D view
%--------------------------------------------------------------------------

%mptPolytope---------------------------------------------------------------
Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2

P1 = polytope(Z1); % convert zonotope Z1 to halfspace representation
P2 = polytope(Z2); % convert zonotope Z2 to halfspace representation

P3 = P1 + P2 % perform Minkowski addition and display result
P4 = P1 & P2; % compute intersection of P1 and P2

V = vertices(P4) % obtain and display vertices of P4

figure
plot(P1); % plot P1
plot(P2); % plot P2
plot(P3,[1 2],'g'); % plot P3
plotFilled(P4,[1 2],[.6 .6 .6],'EdgeColor','none'); % plot P4 
%--------------------------------------------------------------------------

%interval------------------------------------------------------------------
I1 = interval([0; -1], [3; 1]); % create interval I1
I2 = interval([-1; -1.5], [1; -0.5]); % create interval I2
Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1

l = rad(I1) % obtain and display radius of I1
is_intersecting = isIntersecting(I1, Z1) % determines and displays if Z1 is intersecting I1
I3 = I1 & I2; % computes the intersection of I1 and I2
c = mid(I3) % returns and displays the center of I3

figure; hold on
plot(I1); % plot I1
plot(I2); % plot I2
plot(Z1,[1 2],'g'); % plot Z1
plotFilled(I3,[1 2],[.6 .6 .6],'EdgeColor','none'); % plot I3 
%--------------------------------------------------------------------------

%vertices------------------------------------------------------------------
Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
V1 = vertices(Z1); % compute vertices of Z1
A = [0.5 1; 1 0.5]; % numerical matrix A

V2{1} = A*V1; % linear map of vertices
V2{2} = V2{1} + [1; 0]; % translation of vertices
V3 = collect(V2{1},V2); % collect vertices of cell array V2
Zencl = zonotope(V3); % obtain parallelotope containing all vertices

figure
hold on
plot(V2{1},'k+'); % plot V2{1}
plot(V2{2},'ko'); % plot V2{2}
plot(Zencl); % plot Zencl 
%--------------------------------------------------------------------------

%example completed
completed = 1;

%------------- END OF CODE --------------