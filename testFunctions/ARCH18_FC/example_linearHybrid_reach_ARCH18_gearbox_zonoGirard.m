function example_linearHybrid_reach_ARCH18_gearbox_zonoGirard()
% example_linearHybrid_reach_ARCH18_gearbox_zonoGirard -  example of linear
% reachability analysis from the ARCH18 friendly competition (gearbox 
% example). The guard intersections are calculated with Girard's method
%
% Syntax:  
%    example_linearHybrid_reach_ARCH18_gearbox_zonoGirard
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      29-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Options -----------------------------------------------------------------

%initial set
R0 = zonotope([[0; 0; -0.0165; 0.003; 0],diag([0, 0, 1e-4, 1e-4, 0])]);

%initial set
options.x0=center(R0); %initial state for simulation
options.R0=R0; %initial state for reachability analysis

%other
options.startLoc = 1; %initial location
options.finalLoc = 2; %2: gears are meshed
options.tStart=0; %start time
options.tFinal=inf; %final time
options.intermediateOrder = 2;
options.originContained = 0;
options.timeStepLoc{1} = 1.1e-3;
options.timeStepLoc{2} = 1e-3;

options.zonotopeOrder=20; %zonotope order
options.polytopeOrder=3; %polytope order
options.taylorTerms=3;
options.reductionTechnique = 'girard';
options.isHyperplaneMap=0;
options.guardIntersect = 'zonoGirard';
options.enclosureEnables = [1, 2]; %choose enclosure method(s)


%parameters
theta = 0.9; % coefficient of restitution
m_s = 3.2; %[kg], mass of sleve
m_g2 = 18.1; %[kg], mass of second gear
J_g2 = 0.7; %[kg m^2], inertia of second gear
R_s = 0.08; %[m], radius of sleeve
Theta = 36/180*pi; %[rad], included angle of gear
b = 0.01; %[m], width of gear spline
delta_p = -0.002; %[m], p_x sleeve meshes with gear 
%n = 1; % integer number in guard
n = 0; % integer number in guard
F_s = 70; %[N], shifting force
T_f = 1; %[Mm], resisting moment




% Continuous Dynamics -----------------------------------------------------

% system matrix
A = [...
        0 0 0 0 0;...
        0 0 0 0 0;...
        1 0 0 0 0;...
        0 1 0 0 0;...
        0 0 0 0 0]; 
    
B = [1/m_s; 0; 0; 0; 0];

%get dimension
dim = length(A);

%instantiate linear dynamics 
linSys  = linearSys('linearSys',A, eye(dim));




% Hybrid Automaton --------------------------------------------------------

%define large and small distance
dist = 1e3;
smallDist = 1e3*eps;

%first guard set
ch_n = [0, 0, -tan(Theta), -1, 0];
ch_d = -2*n*b;
ch_C(1,:) = [0, 0, 1, 0, 0];
ch_C(2,:) = [0, 0, -1, 0, 0];
ch_C(3,:) = [-sin(Theta), -cos(Theta), 0, 0, 0];
ch_D = [delta_p; b/tan(Theta); 0];
guard1 = constrainedHyperplane(halfspace(ch_n,ch_d),ch_C,ch_D);

C1(1,:) = [0, 0, -tan(Theta), -1, 0];
C1(2,:) = [0, 0, 1, 0, 0];
C1(3,:) = [0, 0, -1, 0, 0];
C1(4,:) = [-sin(Theta), -cos(Theta), 0, 0, 0];
C1(5,:) = -C1(1,:);
d1 = [-2*n*b; delta_p; b/tan(Theta); 0; -2*n*b + smallDist];

%guard = interval([-eps; -dist], [0; -eps]); 
%resets
denominator = m_s*cos(Theta)^2 + m_g2*sin(Theta)^2;
a_11 = (m_s*cos(Theta)^2 - m_g2*theta*sin(Theta)^2)/denominator;
a_12 = (-1*(theta+1)*m_g2*sin(Theta)*cos(Theta))/denominator;
a_21 = (-1*(theta+1)*m_s*sin(Theta)*cos(Theta))/denominator;
a_22 = (m_g2*sin(Theta)^2 - m_s*theta*cos(Theta)^2)/denominator;
a_51 = ((theta+1)*m_s*m_g2*sin(Theta))/denominator;
a_52 = ((theta+1)*m_s*m_g2*cos(Theta))/denominator;
reset1.A = [a_11, a_12, 0, 0, 0; ...
    a_21, a_22, 0, 0, 0; ...]
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    a_51, a_52, 0, 0, 1]; 
reset1.b = zeros(dim,1);

%first transition
trans{1} = transition(guard1,reset1,1,'a','b'); %--> next loc: 1; 'a', 'b' are dummies

%second guard set
ch_n = [0, 0, -tan(Theta), 1, 0];
ch_d = 2*n*b;
ch_C(1,:) = [0, 0, 1, 0, 0];
ch_C(2,:) = [0, 0, -1, 0, 0];
ch_C(3,:) = [-sin(Theta), +cos(Theta), 0, 0, 0];
ch_D = [delta_p; b/tan(Theta); 0];
guard2 = constrainedHyperplane(halfspace(ch_n,ch_d),ch_C,ch_D);

C2(1,:) = [0, 0, -tan(Theta), 1, 0];
C2(2,:) = [0, 0, 1, 0, 0];
C2(3,:) = [0, 0, -1, 0, 0];
C2(4,:) = [-sin(Theta), +cos(Theta), 0, 0, 0];
C2(5,:) = -C2(1,:);
d2 = [2*n*b; delta_p; b/tan(Theta); 0; 2*n*b + smallDist];

%resets
denominator = m_s*cos(Theta)^2 + m_g2*sin(Theta)^2;
a_11 = (m_s*cos(Theta)^2 - m_g2*theta*sin(Theta)^2)/denominator;
a_12 = ((theta+1)*m_g2*sin(Theta)*cos(Theta))/denominator;
a_21 = ((theta+1)*m_s*sin(Theta)*cos(Theta))/denominator;
a_22 = (m_g2*sin(Theta)^2 - m_s*theta*cos(Theta)^2)/denominator;
a_51 = ((theta+1)*m_s*m_g2*sin(Theta))/denominator;
a_52 = (-1*(theta+1)*m_s*m_g2*cos(Theta))/denominator;
reset2.A = [a_11, a_12, 0, 0, 0; ...
    a_21, a_22, 0, 0, 0; ...]
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    a_51, a_52, 0, 0, 1]; 
reset2.b = zeros(dim,1);

%second transition
trans{2} = transition(guard2,reset2,1,'a','b'); %--> next loc: 1; 'a', 'b' are dummies


%third guard set
C3(1,:) = [0, 0, -1, 0, 0];
d3 = [-delta_p]; % different from paper!!
guard3 = halfspace(C3,d3);
%resets
reset3.A = [0, 0, 0, 0, 0; ...
    0, 0, 0, 0, 0; ...
    0, 0, 1, 0, 0; ...
    0, 0, 0, 1, 0; ...
    m_s, m_s, 0, 0, 1]; 
reset3.b = zeros(dim,1);

%third transition
trans{3} = transition(guard3,reset3,2,'a','b'); %--> next loc: 1; 'a', 'b' are dummies

%invariant
Cinv = [-C1(1,:); -C2(1,:); -C1(3,:)];
dinv = [d1(1); d2(1); d1(3)];
dinv = dinv + ones(length(dinv),1)*1e-3;
inv = mptPolytope(Cinv,dinv);
%inv = halfspace([0; 0; 1; 0; 0],-delta_p); %x_3 = p_x 
%inv = interval([-2*eps; -dist], [dist; dist]);

%specify location
loc{1} = location('loc1',1,inv,trans,linSys); 
loc{2} = location('loc2',2,ones(dim,1)*interval(-dist,dist),[],zeroDynSys('zeroDyn',0)); 
%specify hybrid automata
HA = hybridAutomaton(loc); % for "geometric intersection"



% Simulation --------------------------------------------------------------

%set input:
for i=1:2
    options.uLoc{i} = B*F_s + [0; -R_s*T_f/J_g2; 0; 0; 0]; %input for simulation
    options.uLocTrans{i} = options.uLoc{i}; %input center for reachability analysis
    options.Uloc{i} = zonotope(zeros(dim,1)); %input deviation for reachability analysis
end

%simulate hybrid automaton
HA = simulate(HA,options); 




% Reachability Analysis ---------------------------------------------------

%reachable set computations
warning('off','all')
% profile on
tic
[HA] = reach(HA,options);
tComp = toc;
disp(['computation time for gearbox: ',num2str(tComp)]);
% profile viewer
warning('on','all')

disp('Gearbox is meshed')



% Visualization -----------------------------------------------------------

warning('off','all')
figure 
hold on
options.projectedDimensions = [3 4];
options.plotType = 'b';
plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); %plot simulation

% %plot guards
% plot(guard1,[3 4]);
% plot(guard2,[3 4]);
% plot(guard3,[3 4]);

%plot tooth
tooth1 = mptPolytope(C1(1:4,:),d1(1:4));
tooth2 = mptPolytope(C2(1:4,:),d2(1:4));
tooth3 = guard3;
plot(tooth1, [3 4]);
plot(tooth2, [3 4]);
% plot(tooth3, [3 4]);

axis([-2e-2,-2e-3,-1e-2,1e-2]);
box on
warning('on','all')



%------------- END OF CODE --------------