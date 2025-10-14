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

dim_x = 6;

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

params.tStart = 0;
params.tFinal = 7.5;
params.x0 = center(Z_init);
params.R0 = Z_init;
options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 800;
options.reductionTechnique = 'girard';

options.alg = 'lin';
options.tensorOrder = 2;

params.X_sample = zonotope([zeros(dim_x,1),1.8*diag([0.1, 0.1, 1, 0.5, 1, 1])]);


%load reference trajectory
load uTraj_manSpeed

%reference index
refInd = 1597:(1597+750-1); 
u_ref = uTraj(:,refInd);

%obtain parameters and model factors based on friction coefficient mu=1;
p = SRXparameters();
getModelFactors(p,1);

%input consists of reference trajectory u_ref, sensor noise y
params.uTransVec = [u_ref; zeros(4,length(u_ref(1,:)))];
U_ref = zonotope(zeros(length(u_ref(:,1)),1));
Y = zonotope([zeros(4,1), 4*diag([0.02, 0.02, 0.05*pi/180, 0.05*pi/180])]); % applanix data (standard deviation): %delta x,y: 0.02, heading: 0.05*pi/180
params.U = cartProd(U_ref, Y);
options.expVector = [0.2; 0; 0.25; 0.25; 0.15; 0];


% System Dynamics ---------------------------------------------------------

carDyn=nonlinearSys(@DOTBicycleDynamics_controlled_SRX_velEq); 


% Simulation --------------------------------------------------------------

%generate ode options
stepsizeOptions = odeset('MaxStep',params.tFinal-params.tStart);
opt = odeset(stepsizeOptions);

%compute single simulation
inputChanges=ceil(params.tFinal/options.timeStep);
finalTime=params.tFinal;
options.xStep(:,1) = params.x0;

for iChange = 1:inputChanges
    %reset options
    params.tStart=params.tFinal;
    params.tFinal=params.tFinal+finalTime/inputChanges;
    params.u = params.uTransVec(:,iChange);
    if iChange > 1
        params.x0 = x{iChange-1}(:,end);
    end

    %simulate hybrid automaton
    [t{iChange},x{iChange}] = simulate(carDyn,params,opt); 
    options.xStep(:,iChange+1) = x{iChange}(:,end);
end

%reset options
params.tStart = 0;
params.tFinal = finalTime;


% Reachability Analysis ---------------------------------------------------
R = aux_reach(carDyn, params, options);


% Verification ------------------------------------------------------------

%enclose result by interval
IH = interval(R{end});

%saved result
IH_saved = interval( ...
[-0.0476364668386505; 2.6863478436703030; -0.2675705463855188; 64.0741842006564468; 27.8702340052197002; -0.0958282012281481], ...
[-0.0016227388002018; 2.7558057378008112; -0.0343461012859574; 64.4733925792921241; 28.2414003942370044; -0.0057006695447072]);

%final result
assert(isequal(IH,IH_saved,1e-8));


% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function Rcont = aux_reach(sys_init, params, options)

derivatives(sys_init,options); 

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

%obtain intervalhull of uncertain sensor values
IH_u = interval(params.U);
IH_y = IH_u(6:9);

%precompute linearizations
timeSteps = ceil(params.tFinal/options.timeStep);

linParams = cell(1,timeSteps);
linOptions = cell(1,timeSteps);
sys = cell(1,timeSteps);
linsys = cell(1,timeSteps);
RallError = cell(1,timeSteps);

for iSet=1:timeSteps
    
    %set input    
    if iSet>20
        delayInd = 20;
    else
        delayInd = 0;
    end
    params.uTrans = params.uTransVec(:,iSet - delayInd);

    %get first linearized system
    options.linearizationPoint = options.xStep(:,iSet);
    [sysTmp,linsys{iSet},linParams{iSet},linOptions{iSet}] = ...
        linearize(sys_init,[],params,options);
    sys{iSet} = copy(sysTmp);
end

S = struct([]);
for iSet=1:timeSteps
    %prepare reachable set computations
    S = [S; aux_preReach(linsys{iSet},linParams{iSet},linOptions{iSet})];
end
%compute initial additional set due to linearization error
V = zonotope(0*options.expVector,diag(options.expVector));
RallError{1} = particularSolution_timeVarying(linsys{1},...
    V,options.timeStep,options.taylorTerms);

%initialize
Rtp = params.R0;

Rcont = cell(timeSteps,1);
for iSet=1:timeSteps
    
    %set input and center
    if iSet > 30
        delayInd = 30;
    else
        delayInd = 0;
    end
    params.uTrans = params.uTransVec(:,iSet - delayInd);
    c = options.xStep(:,iSet);

    %translate Rinit by linearization point
    Rinit = reduce(Rtp,options.reductionTechnique,options.zonotopeOrder);
    Rdelta = Rinit - c;
    
    % first time step homogeneous solution
    R_hom_tp = S(iSet).eAt*Rdelta + S(iSet).Pu;
    
    %translate reachable sets by linearization point
    R_hom_tp = R_hom_tp + c;

    %do post computations
    [R_hom,IH_hom] = aux_postReach(Rinit,R_hom_tp,...
        S(iSet).PU,S(iSet).F,S(iSet).inputCorr,c);

    %compute maximum reachable set due to maximal allowed linearization error
    IH_max=IH_hom + interval(RallError{iSet});

    % obtain linearization error
    err = linError_constVel(sys{iSet},params.uTrans,IH_max,IH_y);

    %compute performance index of linearization error
    perfInd = max(err./options.expVector);

    % compute reachable set due to the linearization error
    V = zonotope(0*err,diag(err));
    Rerror = particularSolution_timeVarying(linsys{iSet},...
        V,options.timeStep,options.taylorTerms);
    
    %update RallError
    RallError{iSet+1} = enlarge(Rerror,1.8);
    options.expVector = 1.8*err;

    %add intervalhull of actual error
    Rti = R_hom.ti + Rerror;
    Rtp = R_hom.tp + Rerror;
    
    %save reachable set
    Rcont{iSet} = Rti;
end

end

function S = aux_preReach(sys,params,options)
% prepares reachable set computation for linear systems

% compute auxiliary interval matrices
[E,F,G] = taylorMatrices(sys,options.timeStep,options.taylorTerms);

% compute particular solutions
[Pu,inputCorr] = particularSolution_constant(sys,params.uTrans,...
    options.timeStep,options.taylorTerms);
PU = particularSolution_timeVarying(sys,params.U,...
    options.timeStep,options.taylorTerms);

% read out exponential matrix
eAt = getTaylor(sys,'eAdt',struct('timeStep',options.timeStep));

% return in a struct
S = struct('eAt',eAt,'E',E,'F',F,'G',G,'Pu',Pu,'PU',PU,'inputCorr',inputCorr);

end

function [Rnext,IH] = aux_postReach(Rinit,R_tp,RV,F,inputCorr,c)
% computes the reachable continuous set for the first time interval as a
% postprocessing step

%time interval solution
R_err = F*(Rinit+(-c)) + inputCorr;
R_ti = enclose(Rinit,R_tp) + R_err;

%write results to reachable set struct Rfirst
Rnext.tp = R_tp + RV;
Rnext.ti = R_ti + RV;

%compute enclosing hull
IH_init = interval(Rinit);
IH_tp = interval(R_tp);
IH_err = interval(R_err+RV);

IH = or(IH_init,IH_tp) + IH_err;

end

% ------------------------------ END OF CODE ------------------------------
