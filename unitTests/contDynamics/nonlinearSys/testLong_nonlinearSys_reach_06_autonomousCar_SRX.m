function res = testLong_nonlinearSys_reach_06_autonomousCar_SRX()
% testLong_nonlinearSys_reach_06_autonomousCar_SRX - 
%    unit_test_function of nonlinear reachability analysis for following a
%    reference trajectory; the test is similar to
%    test_nonlinear_reach_05_autonomousCar, but here a specialized
%    algorithm is tested that can be run in parallel; the parallelization
%    is not switched on since setting up the threats consumes too much time
%    during unit testing.
%
%    Checks the solution of an autonomous car following a reference
%    trajectory; It is checked whether the reachable set is enclosed in the
%    initial set after a certain amount of time.
%
% Syntax:
%    res = testLong_nonlinearSys_reach_06_autonomousCar_SRX()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       15-March-2012
% Last update:   16-August-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dim_x=6;

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
%initial set
Z_init = zonotope(interval(min(V,[],2),max(V,[],2)));


% Reachability Settings ---------------------------------------------------

options.tStart = 0;
options.tFinal = 7.5;
options.x0 = center(Z_init);
options.R0 = Z_init;
options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 800;
options.reductionTechnique = 'girard';

options.alg = 'lin';
options.tensorOrder = 2;

options.X_sample = zonotope([zeros(dim_x,1),1.8*diag([0.1, 0.1, 1, 0.5, 1, 1])]);


%load reference trajectory
load uTraj_manSpeed

%reference index
refInd = 1597:(1597+750-1); 
u_ref = uTraj(:,refInd);

%obtain parameters and model factors based on friction coefficient mu=1;
p = SRXparameters();
getModelFactors(p,1);

%input consists of reference trajectory u_ref, sensor noise y
options.uTransVec = [u_ref; zeros(4,length(u_ref(1,:)))];
U_ref = zonotope(zeros(length(u_ref(:,1)),1));
Y = zonotope([zeros(4,1), 4*diag([0.02, 0.02, 0.05*pi/180, 0.05*pi/180])]); % applanix data (standard deviation): %delta x,y: 0.02, heading: 0.05*pi/180
options.U = cartProd(U_ref, Y);
options.expVector = [0.2; 0; 0.25; 0.25; 0.15; 0];


% System Dynamics ---------------------------------------------------------

carDyn=nonlinearSys(@DOTBicycleDynamics_controlled_SRX_velEq); 


% Simulation --------------------------------------------------------------

%generate ode options
stepsizeOptions = odeset('MaxStep',options.tFinal-options.tStart);
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
    [t{iChange},x{iChange}] = simulate(carDyn,options,opt); 
    options.xStep(:,iChange+1) = x{iChange}(end,:);
end

%reset options
options.tStart = 0;
options.tFinal = finalTime;


% Reachability Analysis ---------------------------------------------------
Rcont = aux_parallelReach(carDyn, options);


% Verification ------------------------------------------------------------

%enclose result by interval
IH = interval(Rcont{end});

%saved result
IH_saved = interval( ...
[-0.0476364668386505; 2.6863478436703030; -0.2675705463855188; 64.0741842006564468; 27.8702340052197002; -0.0958282012281481], ...
[-0.0016227388002018; 2.7558057378008112; -0.0343461012859574; 64.4733925792921241; 28.2414003942370044; -0.0057006695447072]);

%final result
res = isequal(IH,IH_saved,1e-8);


% Auxiliary functions -----------------------------------------------------

function Rcont = aux_parallelReach(sys_init, options)

derivatives(sys_init,options); 

tic

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

%obtain intervalhull of uncertain sensor values
IH_u = interval(options.U);
IH_y = IH_u(6:9);

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
RallError{1} = errorSolution(linSys{1},linOptions{1},V);

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
    Rinit = reduce(Rtp,options.reductionTechnique,options.zonotopeOrder);
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
    Rerror = errorSolution(linSys{iSet},linOptions{iSet},V);
    
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

% ------------------------------ END OF CODE ------------------------------
