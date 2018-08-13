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
[-5.9905422496220559; -3.3672216670332440; -0.0000000000000011; -0.1358318638060701; -0.2158052418530907; -1.0861985951623208; -111.0558128393433890; -13.4860940275823697; 0.5841357907003362; 1.0310041584072760; -2.9934093674794600; -0.2758782224427850; 0.1112924342645699; -8.5080074621702195; -107.6009707052988773; -2.8508524330906746], ...
[5.5308592400976497; 3.8678774277884425; 0.0000000000000052; 0.2229120840818441; 0.6590247864579211; 0.7065680932148237; 114.7257519415280171; 1.3704587856409969; 0.7150846932613408; 1.2618079267580398; 1.7179670090739867; 0.4806867824268125; 0.6080667597934201; 11.9219058760726035; -7.2684305164878111; 0.7000308978437673]);
        
%check if slightly bloated versions enclose each other
% consider that dim 3 is almost 0
factor = ones(dim,1)*1+1e-8;
factor(3) = 10;
res_1 = (IH <= enlarge(IH_true,factor));
res_2 = (IH_true <= enlarge(IH,factor));

%final result
res = res_1*res_2;

%------------- END OF CODE --------------
