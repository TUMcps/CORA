function example_manual_guard_intersection()
% example_manual_guard_intersection - example from the manual 
% demonstrating the guard intersection
%
% Syntax:
%   example_manual_guard_intersection()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       29-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reconstruction of image, todo real example

% define limits
xLimPt = [15.75,230.8750];
xLimReal = [-0.15,0.20];

yLimPt = [-190.882,-9.382385];
yLimReal = [-4.6,-4.15];

% define guard intersection
guard_V1 = [
 118.5, -162.6703
 149.25, -120.2141
 149.25, -93.26505
 118.5, -93.26505
];
guard_V2 = [
 118.5, -133.5983
 141.0936, -133.5983
 149.25, -121.0121
 149.25, -23.7375
 118.5, -71.18851
];
guard_poly = polygon(guard_V1') | polygon(guard_V2');

% define box enclosure
box_V = [
 118.5, -168.4061 
 149.25, -168.4061
 149.25, -19.20653
 118.5, -19.20653
];
box_poly = polygon(box_V');

% define pca enclosure
pca_V = [
 145.5101, -169.4393
 158.8402, -19.57327
 122.662, -18.1893
 109.3319, -168.0554
];
pca_poly = polygon(pca_V');

% define flow enclosure
flow_V = [
 120.7684, -169.0727
 210.8138, -37.30072
 149.25, -19.20653
 59.20448, -150.9785
];
flow_poly = polygon(flow_V');

% define intersection
intersection_poly = box_poly & pca_poly & flow_poly;

% shift and scale pt to real axis
scale = [xLimReal(2) - xLimReal(1); yLimReal(2) - yLimReal(1)] ...
    ./ [xLimPt(2) - xLimPt(1); yLimPt(2) - yLimPt(1)];
shift = scale .* [sum(xLimPt);sum(yLimPt)]/2 - [sum(xLimReal);sum(yLimReal)]/2;

guard_poly = scale .* guard_poly -shift;
box_poly = scale .* box_poly -shift;
pca_poly = scale .* pca_poly -shift;
flow_poly = scale .* flow_poly -shift;
intersection_poly = scale .* intersection_poly -shift;


% plot --------------------------------------------------------------------

figure;

% plot 1
subplot(1,2,1); hold on; box on;
plot(guard_poly,1:2,'FaceColor',CORAcolor('CORA:reachSet'),'DisplayName','Guard intersection');
plot(box_poly,1:2,'r','DisplayName','box');
plot(pca_poly,1:2,'b','DisplayName','pca');
plot(flow_poly,1:2,'g','DisplayName','flow');
plot(emptySet(2),1:2,'m','DisplayName','Intersection'); % for legend in tikz

xlim(xLimReal); ylim(yLimReal);
xlabel('$x_{(1)}$','Interpreter','latex');
ylabel('$x_{(2)}$','Interpreter','latex');
legend();

% plot 2
subplot(1,2,2); hold on; box on;
plot(guard_poly,1:2,'FaceColor',CORAcolor('CORA:reachSet'));
plot(intersection_poly,1:2,'m');

xlim(xLimReal); ylim(yLimReal);
xlabel('$x_{(1)}$','Interpreter','latex');
% ylabel('$x_{(2)}$','Interpreter','latex');

% ------------------------------ END OF CODE ------------------------------
