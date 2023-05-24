function [HAsim, A, U, U_delay, Uloc, c] = PLL()
% PLL - returns a hybrid model of the phase-locked loop (PLL) described in 
% Sec. 3 of [1]; parameters are from Tab. 1
%
% Syntax:  
%    HAsim = PLL()
%
% Inputs:
%    -
%
% Outputs:
%    HAsim - hybridAutomaton object for simulationg a PLL
%    uLocTrans - cell array of translating inputs
%    U - cell array of input sets
%    U_delay - cell array of input sets for delay phase
%
% References:
%    [1] Althoff, M.; Rajhans, A.; Krogh, B. H.; Yaldiz, S.; Li, X. & Pileggi, L. 
%        Formal Verification of Phase-Locked Loops Using Reachability
%        Analysis and Continuization. Proc. of the Int. Conference on
%        Computer Aided Design, 2011, 659-666.

% Author:       Matthias Althoff
% Written:      10-August-2010
% Last update:  30-June-2020
% Last revision:---

%------------- BEGIN CODE --------------

% parameters
fref=27; %MHz
f0=27e3 - 70; %GHz
N=1000;

Ki=200; %MHz/V 
Kp=25; %MHz/V

%charge pump
Ii=10e-6; %uA
Ip=500e-6; %uA

%capacitors
Ci=25e-12; %pF
Cp1=6.3e-12; %pF
Cp3=2e-12; %pF

%resistors
Rp2=50e3; %kOhm
Rp3=8e3; %kOhm

%states
%x_1 = v_i
%x_2 = v_p1
%x_3 = v_p
%x_4 = phi_v
%x_5 = phi_ref


dim = 6;

%initialize A
A=zeros(dim);

%set values
A(2,2)=-1e-6/Cp1*(1/Rp2+1/Rp3);
A(2,3)=1e-6/(Cp1*Rp3);
A(3,2)=1e-6/(Cp3*Rp3);
A(3,3)=-1e-6/(Cp3*Rp3);
A(4,1)=Ki/N;
A(4,3)=Kp/N;

%obtain input matrix
%number of inputs
nrOfInputs=2;

%initialize B
B=zeros(dim,nrOfInputs);

%set values
B(1,1)=1e-6/Ci;
B(2,2)=1e-6/Cp1;

%obtain constant vector
c = zeros(dim,1);
c(4) = 1/N*f0;
c(5) = fref;
c(6) = 1; %time derivative

% retrieve system dimension
dim_x = length(c);

simSys = linearSys('PLL',A,eye(dim_x)); %initialize system
%--------------------------------------------------------------------------


% create input sets
supFactor = 1.01;
infFactor = 0.99;
vSup = supFactor*[Ii; Ip];
vInf = infFactor*[Ii; Ip];
v = interval(vInf,vSup);
U = B*zonotope(v); 


% create input for delay
diffFactor = supFactor-infFactor;
vSup_delay = diffFactor*[Ii; Ip];
vInf_delay = -diffFactor*[Ii; Ip];
v_delay = interval(vInf_delay,vSup_delay);
U_delay = B*zonotope(v_delay);


% locations:
% loc1: off
% loc2: up
% loc3: dw
% loc4: both 

% define distances
dist = 1e3;
thin = 1e-8;

%loc1:
%input:
Uloc{1} = zonotope(c); 
%invariant
inv = interval([-dist; -dist; -dist; -1-thin; -1-thin; -dist],[dist; dist; dist; 1; 1; dist]);
%guard set 1 sim
guard1Sim = interval([-dist; -dist; -dist; -dist; 1-thin; -dist],[dist; dist; dist; dist; 1; dist]); 
%guard set 2 sim
guard2Sim = interval([-dist; -dist; -dist; 1-thin; -dist; -dist],[dist; dist; dist; 1; dist; dist]);
%reset 1
reset1.A = eye(dim_x);
reset1.c = [0; 0; 0; -1; -1; 0];
%reset 2
reset2.A = eye(dim_x);
reset2.c = [0; 0; 0; -1; -1; 0];
%transition 1 sim
transSim(1) = transition(guard1Sim,reset1,2); %--> next loc: 2
%transition 2 sim
transSim(2) = transition(guard2Sim,reset2,3); %--> next loc: 3
%specify location
locSim(1) = location('off',inv,transSim,simSys);

%loc2:
%input:
Uloc{2} = B*zonotope(v) + c;
%invariant
inv = interval([-dist; -dist; -dist; -1-thin; -1-thin; -dist],[dist; dist; dist; 0; 1; dist]);
%guard set 1 sim
guard1Sim = interval([-dist; -dist; -dist; -thin; -dist; -dist],[dist; dist; dist; 0; dist; dist]); 
%reset 1
Atmp = eye(dim_x);
Atmp(6,6) = 0;
reset1.A = Atmp;
reset1.c = zeros(dim_x,1);
%transition 1 sim
clear transSim
transSim = transition(guard1Sim,reset1,4); %--> next loc: 4
%specify location
locSim(2) = location('up',inv,transSim,simSys);


%loc3:
%input:
Uloc{3} = -B*zonotope(v) + c;
%invariant
inv = interval([-dist; -dist; -dist; -1-thin; -1-thin; -dist],[dist; dist; dist; 1; 0; dist]);
%guard set 1 sim
guard1Sim = interval([-dist; -dist; -dist; -dist; -thin; -dist],[dist; dist; dist; dist; 0; dist]);
%reset 1
Atmp = eye(dim_x);
Atmp(6,6) = 0;
reset1.A = Atmp;
reset1.c = zeros(dim_x,1);
%transition 1 sim
transSim = transition(guard1Sim,reset1,4); %--> next loc: 4
%specify location
locSim(3) = location('dw',inv,transSim,simSys);


%loc4:
%input:
Uloc{4} = B*zonotope(v_delay) + c;
%invariant
inv = interval([-dist; -dist; -dist; -dist; -dist; -dist],[dist; dist; dist; dist; dist; dist]);
%guard set 1 sim
guard1Sim = interval([-dist; -dist; -dist; -dist; -dist; 50e-6],[dist; dist; dist; dist; dist; 51e-6]);
%reset 1
reset1.A = eye(dim_x);
reset1.c = zeros(dim_x,1);
%transition 1 sim
transSim = transition(guard1Sim,reset1,1); %--> next loc: 1
%specify location
locSim(4) = location('both',inv,transSim,simSys);

% instantiate hybrid automaton
HAsim = hybridAutomaton(locSim);

