function res = test_nonlinear_reach_06_autonomousCarOptimized_SRX()
% test_nonlinear_reach_06_autonomousCarOptimized_SRX - unit_test_function of 
% nonlinear reachability analysis for following a reference trajectory; the
% test is similar to test_nonlinear_reach_05_autonomousCar, but here a
% specialized algorithm is tested that can be run in parallel; the
% parallelization is not switched on since setting up the threats consumes
% to much time during unit testing.
%
% Checks the solution of an autonomous car following a reference
% trajectory;
% It is checked whether the reachable set is enclosed in the initial set
% after a certain amount of time.
%
% Syntax:  
%    res = test_nonlinear_reach_06_autonomousCarOptimized_SRX()
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
% Written:      15-March-2012
% Last update:  16-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

dim=6;

%load data
load SRX_lc_data

%obtain interval hull of initial states
for i=1:13
    V(1,i) = 0; % slip
    V(2,i) = x_cad{i}(2,1); % yaw angle
    V(3,i) = x_cad{i}(3,1); % yaw rate
    V(4,i) = x_cad{i}(4,1); % x-pos
    V(5,i) = x_cad{i}(5,1); % y-pos
    V(6,i) = x_cad{i}(6,1); % wheel angle
end
vert = vertices(V);
%initial set
Z_init = zonotope(interval(vert));

%set options --------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=7.5; %start time
%options.tFinal=3; %start time
options.x0=center(Z_init); %initial state for simulation
options.R0=Z_init; %initial state for reachability analysis
options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
%options.zonotopeOrder=120; %zonotope order
options.zonotopeOrder=800; %zonotope order
options.polytopeOrder=1; %polytope order

options.originContained = 0;
options.tensorOrder = 1;
options.intermediateOrder=[];

options.X_sample = zonotope([zeros(dim,1),1.8*diag([0.1, 0.1, 1, 0.5, 1, 1])]);
%--------------------------------------------------------------------------

%load reference trajectory
load uTraj_manSpeed

%reference index
refInd = 1597:(1597+750-1); 
u_ref = uTraj(:,refInd);

%obtain parameters and model factors based on friction coefficient mu=1;
p = SRXparameters();
getModelFactors(p,1);

%input consists of reference trajectory u_ref, sensor noise y
options.uTransVec = [u_ref; zeros(5,length(u_ref(1,:)))];
U_ref = zonotope(zeros(length(u_ref(:,1)),1));
Y = zonotope([zeros(5,1), 4*diag([0.02, 0.02, 0.05*pi/180, 0.05*pi/180, 0.02])]); % applanix data (standard deviation): %delta x,y: 0.02, heading: 0.05*pi/180
options.U = cartesianProduct(U_ref, Y);
options.expVector = [0.2; 0; 0.25; 0.25; 0.15; 0];

%specify continuous dynamics-----------------------------------------------
carDyn=nonlinearSys(6,10,@DOTBicycleDynamics_controlled_SRX_velEq,options); %initialize car
%--------------------------------------------------------------------------


%generate ode options
stepsizeOptions = odeset('MaxStep',options.tStart-options.tFinal);
opt = odeset(stepsizeOptions);

%compute single simulation
inputChanges=ceil(options.tFinal/options.timeStep);
finalTime=options.tFinal;
options.xStep(:,1) = options.x0;

for iChange = 1:inputChanges
    %reset options
    options.tStart=options.tFinal;
    options.tFinal=options.tFinal+finalTime/inputChanges;
    options.u = options.uTransVec(:,iChange);
    if iChange > 1
        options.x0 = x{iChange-1}(end,:);
    end

    %simulate hybrid automaton
    [carDyn,t{iChange},x{iChange}] = simulate(carDyn,options,options.tStart,options.tFinal,options.x0,opt); 
    options.xStep(:,iChange+1) = x{iChange}(end,:);
end

%reset options
options.tStart = 0;
options.tFinal = finalTime;

%reachable set for no uncertain parameters
Rcont = parallelReach(carDyn, options);

%enclose result by interval
IH = interval(Rcont{end});

%saved result
IH_saved = interval( ...
    [-0.0476364710979309; 2.6863478250960333; -0.2675705811035878; 64.0741841510284615; 27.8702338802093976; -0.095828222990821], ...
    [-0.0016227345409213; 2.7558057563750809; -0.0343460665678876; 64.4733926289201094; 28.2414005192472999; -0.005700647782034]);

%check if slightly bloated versions enclose each other
res_1 = (IH <= enlarge(IH_saved,1+1e-8));
res_2 = (IH_saved <= enlarge(IH,1+1e-8));

%final result
res = res_1*res_2;



function Rcont = parallelReach(sys_init, options)


tic

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end

%obtain intervalhull of uncertain sensor values
IH_u = interval(options.U);
IH_y = IH_u(6:10);

%precompute linearizations
timeSteps = ceil(options.tFinal/options.timeStep);

linOptions=cell(1,timeSteps);
sys=cell(1,timeSteps);
linSys=cell(1,timeSteps);
RallError=cell(1,timeSteps);


for iSet=1:timeSteps
    
    %set input    
    if iSet>20
        delayInd = 20;
    else
        delayInd = 0;
    end
    options.uTrans = options.uTransVec(:,iSet - delayInd);

    %get first linearized system
    options.linearizationPoint = options.xStep(:,iSet);
    [sysTmp,linSys{iSet},linOptions{iSet}] = linearize(sys_init,options);
    sys{iSet} = copy(sysTmp);
end
    
%parfor iSet=1:timeSteps
for iSet=1:timeSteps  
    %prepare reachable set computations
    linSys{iSet} = preReach(linSys{iSet}, linOptions{iSet});
end
%compute initial additional set due to linearization error
V = zonotope([0*options.expVector,diag(options.expVector)]);
RallError{1} = errorSolution(linSys{1},V,linOptions{1});

%initialize
t = options.tStart + options.timeStep;
iSet = 1;
perfInd = 0;
Rtp = options.R0;

toc

tic

%while final time is not reached
while (t<options.tFinal) && perfInd<1
    
    %set input and center
    if iSet>30
        delayInd = 30;
    else
        delayInd = 0;
    end
    options.uTrans = options.uTransVec(:,iSet - delayInd);
    options.center = options.xStep(:,iSet);

    %translate Rinit by linearization point
    Rinit = reduce(Rtp,'girard',options.zonotopeOrder);
    Rdelta = Rinit + (-options.center);

    %do core computations
    R_hom_tp = coreReach(linSys{iSet}, Rdelta);
    
    %translate reachable sets by linearization point
    R_hom_tp = R_hom_tp+options.center;

    %do post computations
    [R_hom,IH_hom] = postReach(linSys{iSet}, Rinit, R_hom_tp, options.center);

    %compute maximum reachable set due to maximal allowed linearization error
    IH_max=IH_hom + interval(RallError{iSet});

    % obtain linearization error
    [error] = linError_constVel(sys{iSet},options.uTrans,IH_max,IH_y);

    %compute performance index of linearization error
    perfInd = max(error./options.expVector)
    
    if perfInd>1
        disp('stop');
    end

    % compute reachable set due to the linearization error
    V = zonotope([0*error,diag(error)]);
    [Rerror] = errorSolution(linSys{iSet},V,linOptions{iSet});
    
    %update RallError
    RallError{iSet+1} = enlarge(Rerror,1.8);
    options.expVector = 1.8*error;

    %add intervalhull of actual error
    Rti =R_hom.ti + Rerror;
    Rtp =R_hom.tp + Rerror;
    
    %save reachable set
    Rcont{iSet}=Rti; 
    
    %increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1; 
    %t
    
end
toc

%------------- END OF CODE --------------
