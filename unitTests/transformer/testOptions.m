% Text = 273.15;
Text = 0;

% options
%set options---------------------------------------------------------------
options.x0 = [Text+16; Text+16]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([1, 1])]); %initial set
options.startLoc = [1 1]; %initial location
options.finalLoc = 0; %0: no final location
options.tStart = 0; %start time
options.tFinal = 24 * 3600; %final time
options.timeStepLoc{1} = 1; %time step size in location 1
options.timeStepLoc{2} = 1; %time step size in location 1
options.timeStepLoc{3} = 1; %time step size in location 1
options.timeStepLoc{4} = 1; %time step size in location 1
options.taylorTerms = 10; % reachability ?
options.polytopeType = 'mpt';

 options.zonotopeOrder = 20;
 options.polytopeOrder = 10;
 options.errorOrder=2;
 options.reductionTechnique = 'girard';
 options.isHybrid = 1;
 options.isHyperplaneMap = 0;
 options.enclose = [5]; %choose enclosure method(s)
 options.originContained = 0;
%--------------------------------------------------------------------------
%set input:
options.uLoc{1} = [1; 1]; %input for simulation
options.uLocTrans{1} = options.uLoc{1}; %center of input set
options.Uloc{1} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{2} = [1; 1]; %input for simulation
options.uLocTrans{2} = options.uLoc{2}; %center of input set
options.Uloc{2} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{3} = [1; 1]; %input for simulation
options.uLocTrans{3} = options.uLoc{3}; %center of input set
options.Uloc{3} = zonotope(zeros(2,1)); %input deviation from center

options.uLoc{4} = [1; 1]; %input for simulation
options.uLocTrans{4} = options.uLoc{4}; %center of input set
options.Uloc{4} = zonotope(zeros(2,1)); %input deviation from center