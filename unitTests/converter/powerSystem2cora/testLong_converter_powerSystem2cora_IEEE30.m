function res = testLong_converter_powerSystem2cora_IEEE30()
% testLong_converter_powerSystem2cora_IEEE30 - unit test for 
%    converting an IEEE 30 bus benchmark problem in [1] to a CORA model
%    It is checked whether the dynamic and constraint functions match for 
%    random variables
%
% Syntax:
%    res = testLong_converter_powerSystem2cora_IEEE30()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Authors:       Matthias Althoff
% Written:       17-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% create CORA models
% full system
powerSystem2cora('IEEE30')
% subsystem 1
powerSystem2cora('IEEE30_sub1')
% subsystem 2
powerSystem2cora('IEEE30_sub2')

%% load converted models
% set path
path = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];
% load full model
load([path filesep 'IEEE30_model'], 'IEEE30_model');
% load subsystem 1
load([path filesep 'IEEE30_sub1_model'], 'IEEE30_sub1_model');
% load subsystem 2
load([path filesep 'IEEE30_sub2_model'], 'IEEE30_sub2_model');

% rewrite systems as cells
sys{1} = IEEE30_model;
sys{2} = IEEE30_sub1_model;
sys{3} = IEEE30_sub2_model;

%% create saved dynamics
% create full model
IEEE30_groundTruth_model = nonlinDASys('IEEE30_groundTruth',@IEEE30_groundTruth_dyn,@IEEE30_groundTruth_con,18,36,59);
% create subsystem 1
IEEE30_groundTruth_sub1_model = nonlinDASys('IEEE30_groundTruth_sub1',@IEEE30_groundTruth_sub1_dyn,@IEEE30_groundTruth_sub1_con,9,30,33);
% create subsystem 2
IEEE30_groundTruth_sub2_model = nonlinDASys('IEEE30_groundTruth_sub2',@IEEE30_groundTruth_sub2_dyn,@IEEE30_groundTruth_sub2_con,9,28,26);

% rewrite systems as cells
sys_saved{1} = IEEE30_groundTruth_model;
sys_saved{2} = IEEE30_groundTruth_sub1_model;
sys_saved{3} = IEEE30_groundTruth_sub2_model;

% inputs of the full model
u_gen = [2.6; 0.4; 0; 0; 0; 0]; %from U Washington website
u_load = zeros(sys{1}.nrOfInputs - length(u_gen),1);
u0full = [u_gen; u_load];

%% nominal values for the input, state, and constraints of the full system
% obtain number of generators and number of buses
nrOfGenerators = length(u_gen);
nrOfBuses = length(u_load);

% reasonable initial state (dynamic and algebraic)
x0full = [0*ones(nrOfGenerators,1); 377*ones(nrOfGenerators,1); 1*ones(nrOfGenerators,1)];
y0full = [ones(nrOfBuses, 1); zeros(nrOfBuses-1, 1)];

% derivatives to update jacobian
options.tensorOrder = 2;
derivatives(IEEE30_model,options);

% compute initial state solution of the system and the saved system
y0full = consistentInitialState(sys{1}, x0full, y0full, u0full);

%% build indices for the input, state, and constraints relating the full 
%% system to the subsystems
% full system
load IEEE30.mat bus
busFull = bus;

% subsystem 1
load IEEE30_sub1.mat bus
bus_subsystem{1} = bus;

% subsystem 2
load IEEE30_sub2.mat bus
bus_subsystem{2} = bus;

% create cell of subsystems
subsystem{1} = sys{2};
subsystem{2} = sys{3};

% compute indices for subsystems
index = indexForSubsystems(sys{1}, subsystem, busFull, bus_subsystem, 3);

% set accuracy
accuracy = 1e-6;

% initialize partial results
resPartial = [];


%% perform evaluation
% loop over systems
for i = 1:length(sys)   
    % use indices to map states and inputs of the full system
    if i>1 % subsystem nr starts at 2
        % project initial dynamic state
        x0 = index{i-1}.X*x0full;
        % project initial algebraic state
        y0 = index{i-1}.Y*y0full;
        % project inputs
        u0 = index{i-1}.Uu*u0full;
    else
        x0 = x0full;
        y0 = y0full;
        u0 = u0full;
    end

    for iRun = 1:10
        % create random values for deviations of states, inputs, and constraints
        x = x0 + (0.2*rand(sys{i}.dim,1) - 0.1*ones(sys{i}.dim,1)).*x0; %+-10% of steady state
        y = y0 + (0.2*rand(sys{i}.nrOfConstraints,1) - 0.1*ones(sys{i}.nrOfConstraints,1)).*y0; %+-10% of steady state
        u = u0 + (0.2*rand(sys{i}.nrOfInputs,1) - 0.1*ones(sys{i}.nrOfInputs,1)).*u0; %+-10% of steady state
        
        % dynamic function output
        f = sys{i}.dynFile(x,y,u);
        f_saved = sys_saved{i}.dynFile(x,y,u);
        
        % constraint function output
        g = sys{i}.conFile(x,y,u);
        g_saved = sys_saved{i}.conFile(x,y,u);
        
        %check if slightly bloated versions enclose each other
        resPartial(end+1) = (norm(f - f_saved) <= accuracy);
        resPartial(end+1) = (norm(g - g_saved) <= accuracy);
    end
end

% final result
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
