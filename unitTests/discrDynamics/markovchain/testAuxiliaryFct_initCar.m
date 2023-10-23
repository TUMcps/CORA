function [HA,options,params,stateField,inputField,changeSpeed] = testAuxiliaryFct_initCar(varargin)
% testAuxiliaryFct_initCar - initializes a car model for the abstraction to
%    a Markov chain
%
% Syntax:
%    [HA,options,stateField,inputField,changeSpeed] = testAuxiliaryFct_initCar()
%
% Inputs:
%    -
%
% Outputs:
%    HA - hybrid automaton object of the car
%    options - options structure for the simulation of the car
%    stateField - partition object of the state space
%    inputField - partition object of the input space
%    changeSpeed - speed from which the acceleration is bounded by engine
%    power

% Authors:       Matthias Althoff
% Written:       12-October-2009
% Last update:   31-July-2016
%                09-July-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% set parameters
params.tStart = 0; %start time
params.tFinal = 0.5; %final time
params.startLoc = 1; %initial location

%% set options
options.timeStep = 0.05; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.zonotopeOrder = 10; %zonotope order
options.polytopeOrder = 6; %polytope order
options.target = 'vehicleDynamics';

%% specify continuous dynamics
accSlow = linearSys('accEidSlow',[0 1;0 0],[0;7]); %acceleration
accFast = nonlinearSys('accEidFast',@accSysEidFast); %acceleration
dec = linearSys('decEid',[0 1;0 0],[0;7]); %deceleration
sL = linearSys('sL',[0 1;0 0],[0;0]); %speed limit
sS = linearSys('sS',[0 0;0 0],[0;0]); %standstill


%% specify transitions

%reset map for all transitions
reset.A = eye(2); 
reset.b = zeros(2,1);

%values required to set up invariant and guard sets
dist = 1e3; 
eps = 1e-6;
maxSpeed = 20;
changeSpeed = 7.3;

%specify invariant
inv = interval([0; 0], [dist; maxSpeed]); %invariant for all locations

%guard sets
Istop = interval([0; 0], [dist; eps]);
ImaxSpeed = interval([0; maxSpeed-eps], [dist; maxSpeed]);  
IaccChange = interval([0; changeSpeed], [dist; changeSpeed+eps]); 

%specify transitions
tran1{1} = transition(IaccChange,reset,2); 
tran2{1} = transition(ImaxSpeed,reset,3); 
tran3 = []; 
tran4{1} = transition(Istop,reset,5);
tran5 = [];

%% specify locations              
loc{1} = location('accSlow',inv,tran1,accSlow);
loc{2} = location('accFast',inv,tran2,accFast);
loc{3} = location('sL',inv,tran3,sL);
loc{4} = location('dec',inv,tran4,dec);
loc{5} = location('sS',inv,tran5,sS);


%% specify hybrid automaton
HA = hybridAutomaton(loc);

%% initialize partition
stateField = partition([0, 200;... %position in m
                      0, 20],... %velocity in m/s         
                     [40;10]);
inputField = partition([-1,1],...  %acceleartion in m/s^2
                     6);  
end

% ------------------------------ END OF CODE ------------------------------
