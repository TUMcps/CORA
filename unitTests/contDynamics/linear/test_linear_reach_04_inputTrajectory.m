function res = test_linear_reach_04_inputTrajectory()
% test_linear_reach_04_inputTrajectory - unit test for linear reachability 
% analysis with an input trajectory uTransVec; this test should check 
% whether correct the input trajectory is correctly considered
%
% Syntax:  
%    res = test_linear_reach_04_inputTrajectory()
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
% Written:      27-July-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% load data
load data_Jean-MarcBiannic.mat A B uvec
dim = length(A);
inputDim = length(B(1,:));

%set options --------------------------------------------------------------
options.tStart = 0; %start time
options.tFinal = 10; %final time
options.x0 = zeros(dim,1); %initial state for simulation
options.R0 = zonotope([options.x0,0.1*eye(dim,4)]); %initial state for reachability analysis

options.timeStep = 0.01; %time step size for reachable set computation
options.taylorTerms = 4; %number of taylor terms for reachable sets
options.zonotopeOrder = 50; %zonotope order
options.originContained = 0;
options.reductionTechnique = 'girard';
options.linAlg = 1;

options.uTransVec=[zeros(inputDim,1) uvec]; % input trajectory
options.U = zonotope([zeros(inputDim,1),diag([0.05 1])]); %input for reachability analysis

options.linAlg=2;
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
sys = linearSys('JeanMarcSys',A,B); %initialize system
%--------------------------------------------------------------------------

%compute reachable set using zonotopes
Rcont = reach(sys, options);

IH = interval(Rcont{end});

%saved result
IH_true = interval( ...
[-5.9905560293201070; -3.4057356893286905; -0.0000000000000005; -0.1358299246418150; -0.2158580931624483; -1.0965659219321691; -112.1115569589848633; -13.5901685641709289; 0.5841279576971828; 1.0309215393525626; -3.0059985516073930; -0.2762748408328315; 0.1112566862281668; -8.5076833881716656; -107.5987924139651852; -2.8870004572813155], ...
[5.5371951163933488; 3.8685515060073858; 0.0000000000000050; 0.2242871135888930; 0.6633474461607697; 0.7066935454436742; 117.4322678652501111; 1.4059728502365960; 0.7150847122292477; 1.2618079325715417; 1.7199593172943310; 0.4827850307499247; 0.6080613032843966; 11.9394930370437056; -7.0395105967679683; 0.7394527859620696]);
        
%check if slightly bloated versions enclose each other
% consider that dim 3 is almost 0
factor = ones(dim,1)*1+1e-8;
factor(3) = 10;
res_1 = (IH <= enlarge(IH_true,factor));
res_2 = (IH_true <= enlarge(IH,factor));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
