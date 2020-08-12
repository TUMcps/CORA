function [HA,options,stateField,inputField]=initBicycle(varargin)
% changed: 02-November-2009


directory=cd;

%set options---------------------------------------------------------------
options.tStart=0; %start time
options.tFinal=0.5; %final time
%options.tFinal=0.5; %final time
options.startLoc=1; %initial location
options.finalLoc=0; %final location
options.timeStep=0.05; %time step size for reachable set computation
options.taylorTerms=4; %number of taylor terms for reachable sets
options.zonotopeOrder=10; %zonotope order
options.polytopeOrder=6; %polytope order
options.projectedDimensions=[1 2];
options.reductionInterval=1e3;
options.target='vehicleDynamics';
options.guardIntersect='polytope';
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
accSlow=linearSys('accEidSlow',[0 1;0 0],[0;2.5]); %acceleration
accFast=nonlinearSys('accEidBicycle',@accSysEidFastBicycle); %acceleration
dec=linearSys('decEid',[0 1;0 0],[0;7]); %deceleration
sL=linearSys('sL',[0 1;0 0],[0;0]); %speed limit
sS=linearSys('sS',[0 0;0 0],[0;0]); %standstill
%--------------------------------------------------------------------------

%specify transitions-------------------------------------------------------

%reset map for all transitions
reset.A=eye(2); 
reset.b=zeros(2,1);

%long distance l
l=[0 1e3]; %-->probabilities may only escape in driving direction
maxSpeed=[19 19.99];
zeroSpeed=[0.001 0.01];
accChangeSpeed=[2.5 2.55];

%specify invariant
inv=interval([l; 0 20]); %invariant for all locations

%guard sets
IHstop=interval([l; zeroSpeed]);
IHmaxSpeed=interval([l; maxSpeed]);  
IHaccChange=interval([l; accChangeSpeed]); 

%specify transitions
tran1{1}=transition(IHaccChange,reset,2); 
tran2{1}=transition(IHmaxSpeed,reset,3); 
tran3=[]; 
tran4{1}=transition(IHstop,reset,5);
tran5=[];

%--------------------------------------------------------------------------

%specify locations              
loc{1}=location('accSlow',inv,tran1,accSlow);
loc{2}=location('accFast',inv,tran2,accFast);
loc{3}=location('sL',inv,tran3,sL);
loc{4}=location('dec',inv,tran4,dec);
loc{5}=location('sS',inv,tran5,sS);


%specify hybrid automaton
HA=hybridAutomaton(loc);

%Initialize partition------------------------------------------------------
stateField=partition([0, 200;... %position in m
                      0, 20],... %velocity in m/s         
                     [40;10]);
inputField=partition([-1,1],...  %acceleartion in m/s^2
                     6);  
%--------------------------------------------------------------------------

if nargin==1
    cd(directory);
end