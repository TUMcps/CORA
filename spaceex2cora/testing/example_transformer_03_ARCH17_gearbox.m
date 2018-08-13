function example_transformer_03_ARCH17_gearbox()
% example_linearHybrid_reach_ARCH17_gearbox -  example of linear reachability 
% analysis from the ARCH17 friendly competition (gearbox example); the
% linear dynamics switches according to guard sets
%
% Syntax:  
%    example_linearHybrid_reach_ARCH17_gearbox
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
% Author:       Matthias Althoff
% Written:      17-March-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
options.enclosureEnables = [3, 5]; %choose enclosure method(s)

%build options for SX automaton
options_SX = options;

%add extra dimension [t]
R0_SX = zonotope([[0; 0; 0; -0.0165; 0.003; 0],diag([0, 0, 0, 1e-4, 1e-4, 0])]);
options_SX.x0 = center(R0_SX);
options_SX.R0 = R0_SX;

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
% Continuous dynamics------------------------------------------------------
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
%--------------------------------------------------------------------------

%specify hybrid automaton--------------------------------------------------
%define large and small distance
dist = 1e3;
smallDist = 1e3*eps;

%first guard set
C1(1,:) = [0, 0, -tan(Theta), -1, 0];
C1(2,:) = [0, 0, 1, 0, 0];
C1(3,:) = [0, 0, -1, 0, 0];
C1(4,:) = [-sin(Theta), -cos(Theta), 0, 0, 0];
C1(5,:) = -C1(1,:);
d1 = [-2*n*b; delta_p; b/tan(Theta); 0; -2*n*b + smallDist];
guard1 = mptPolytope(C1,d1);
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
C2(1,:) = [0, 0, -tan(Theta), 1, 0];
C2(2,:) = [0, 0, 1, 0, 0];
C2(3,:) = [0, 0, -1, 0, 0];
C2(4,:) = [-sin(Theta), +cos(Theta), 0, 0, 0];
C2(5,:) = -C2(1,:);
d2 = [2*n*b; delta_p; b/tan(Theta); 0; 2*n*b + smallDist];
guard2 = mptPolytope(C2,d2);
%guard = interval([-eps; -dist], [0; -eps]); 
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
guard3 = mptPolytope(C3,d3);
%guard = interval([-eps; -dist], [0; -eps]); 
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

%build hybrid automaton from SX file
%sha_gears = SX2structHA('SX_Mesh.xml','mesh');
%StructHA2file(sha_gears,'gearboxSX');
HA_SX = gearboxSX();
%--------------------------------------------------------------------------

%set input:
for i=1:2
    options.uLoc{i} = B*F_s + [0; -R_s*T_f/J_g2; 0; 0; 0]; %input for simulation
    options.uLocTrans{i} = options.uLoc{i}; %input center for reachability analysis
    options.Uloc{i} = zonotope(zeros(dim,1)); %input deviation for reachability analysis
    
    %SX model only has 1 dummy input (has no effect)
    options_SX.uLoc{i} = 0; %input for simulation
    options_SX.uLocTrans{i} = options_SX.uLoc{i}; %input center for reachability analysis
    options_SX.Uloc{i} = zonotope(0); %input deviation for reachability analysis
end

%simulate hybrid automata
HA = simulate(HA,options);
temp = get(HA,'result');
simRes = temp.simulation;

HA_SX = simulate(HA_SX,options_SX); 
temp = get(HA_SX,'result');
simRes_SX = temp.simulation;

%reachable set computations
%profile on
% REACHABILITY OF HA_SX SEEMS TO NOT TERMINATE
% tic
% [HA] = reach(HA,options);
% tComp = toc;
% disp(['computation time for gearbox: ',num2str(tComp)]);
% 
% tic
% [HA_SX] = reach(HA_SX,options_SX);
% tComp = toc;
% disp(['computation time for gearbox: ',num2str(tComp)]);
%profile viewer

disp('Gearbox is meshed')

% COMPUTE UNIT TEST SUCCESS
% check number of timesteps
t1 = simRes.t{1};
t2 = simRes_SX.t{1};
tSteps = min(length(t1),length(t2));

% evaluate difference of variables px & py
px1 = simRes.x{1}(:,3);
px2 = simRes_SX.x{1}(:,4);
py1 = simRes.x{1}(:,4);
py2 = simRes_SX.x{1}(:,5);

nx = norm(px1(1:tSteps) - px2(1:tSteps));
ny = norm(py1(1:tSteps) - py2(1:tSteps));
if(nx + ny > 1e-1)
    error('unit test failed, error=(%i,%i)',nx,ny);
else
    disp("error=(" + string(nx) + "," + string(ny) + "), good job!");
end


%choose projection and plot------------------------------------------------
figure 
hold on
options.projectedDimensions = [3 4];
options.plotType = 'b';
% plot(HA,'reachableSet',options); %plot reachable set
plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA,'simulation',options); %plot simulation

options_SX.projectedDimensions = [4 5];
options_SX.plotType = 'b';
plotFilled(options_SX.R0,options_SX.projectedDimensions,'w','EdgeColor','k'); %plot initial set
plot(HA_SX,'simulation',options_SX); %plot simulation

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
plot(tooth3, [3 4]);

axis([-2e-2,-2e-3,-1e-2,1e-2]);
%--------------------------------------------------------------------------


%------------- END OF CODE --------------