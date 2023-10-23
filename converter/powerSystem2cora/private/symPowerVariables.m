function [powVariables, NrOf] = symPowerVariables(scenario)
% symPowerVariables - create the following symbolic variables: 
% 1 is the slack bus
% 2...(nrOfGenBuses+1) are the generator buses
% (nrOfGenBuses+2)...(nrOfGenBuses+nrOfLoadBuses+2) are the load buses
%
% dynamic states:
% x(i) = delta(i), i = 1...nrOfGenerators; generator phase angles
% x(nrOfGenerators + i) = omega(i), i = 1...nrOfGenerators; angular velocity
% x(2*nrOfGenerators + i) = P_m(i), i = 1...nrOfGenerators; mechanical
% power
%
% algebraic states:
% y(i) = E(i), i = 1...nrOfGenerators; generator voltage
% y(nrOfGenerators + i) = V(nrOfGenerators + i), i = 1...nrOfLoadBuses; bus voltage
% y(nrOfGenerators + nrOfLoadBuses + i - 1) = Theta(i) - delta, j = 2...nrOfBuses (first bus is the slack bus)
%
% inputs:
% u(i) = P_c(i), i = 1...nrOfGenerators; commanded power production
% u(nrOfGenerators + i) = P^d_g(i), i = 1...nrOfBuses; directly injected active power
% u(nrOfGenerators + nrOfBuses + i) = Q^d_g(i), i = 1...nrOfBuses; directly injected reactive power
% u(nrOfGenerators + 2*nrOfBuses + i) = V(i), i = 1...nrOfInputs; input voltages
% u(nrOfGenerators + 2*nrOfBuses + nrOfInputs + i) = Theta(i), i = 1...nrOfInputs; input phases
%
% More information can be found in [1, Sec. VII].
%
% Syntax:
%    [powVariables, NrOf] = symPowerVariables(scenario)
%
% Inputs:
%    scenario - struct specifying a power system scenario
%
% Outputs:
%    powVariables - symbolic variables for the power system
%    NrOf - ???
%
% References:
%    [1] M. Althoff, "Benchmarks for the Formal Verification of Power 
%        Systems", Proc. of the 9th International Workshop on Applied 
%        Verification of Continuous and Hybrid Systems, 
%        2022, x-x

% Authors:       Matthias Althoff
% Written:       14-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain bus struct
bus = scenario.bus;

%% number of buses (NoB)
NoB.input = length(bus.input); % nr. of buses serving as inputs
NoB.generator = length(bus.generator); % nr. of generators
NoB.load = length(bus.load); % nr. of loads
NoB.total = length(scenario.Y) - NoB.input; % total number of buses

% sanity check
if NoB.total~=(NoB.generator + NoB.load)
    disp('Number of buses is inconsistent');
end

%% number of states, inputs, and constraints
NrOf.states = 3*NoB.generator; % each generator bus has three state variables
NrOf.inputs = NoB.generator + 2*NoB.input + NoB.total;
if ~isempty(bus.slack) % in case of slack bus, remove one constraint variable
    NrOf.constraints = 2*NoB.total - 1;
else
    NrOf.constraints = 2*NoB.total;
end

% generate symbolic constraint states
for i=1 : NrOf.constraints
    command=['y(',num2str(i),',1)=sym(''yL',num2str(i),'R'');'];
    eval(command);
end 

% generate symbolic dynamic states
for i=1 : NrOf.states
    command=['x(',num2str(i),',1)=sym(''xL',num2str(i),'R'');'];
    eval(command);
end 

% generate inputs
for i=1 : NrOf.inputs
    command=['u(',num2str(i),',1)=sym(''uL',num2str(i),'R'');'];
    eval(command);
end 

%% set generator voltages
E = y(1:NoB.generator);
% faulty bus exists
if ~isempty(bus.fault)
    % set voltage at faulty bus to 0
    E(bus.fault, 1) = sym(0);
end

%% set voltages for slack and generator buses
% generator buses exist
if NoB.generator>0
    V(1:NoB.generator, 1) = sym(scenario.VM); 
    % faulty bus exists
    if ~isempty(bus.fault) 
        % set voltage to faulty bus to E(1)
        V(bus.fault, 1) =  y(1); 
    end
end

%% set voltages for load buses
V(NoB.generator + 1 : NoB.total, 1) = y(NoB.generator + 1 : NoB.total, 1); %internal voltages
V(NoB.total + 1 : NoB.total + NoB.input, 1) = u(NoB.generator + 1 : NoB.generator + NoB.input, 1); %external voltages

%% set phase angles
% subsystem contains slack bus
if ~isempty(bus.slack)
    Theta(1, 1) = sym(0);
    Theta(2:NoB.total, 1) = y(NoB.total + 1 : 2*NoB.total - 1, 1); %internal phases
% subsystem does not contain the slack bus
else
    Theta(1:NoB.total, 1) = y(NoB.total + 1 : 2*NoB.total, 1); %internal phases
end
Theta(NoB.total + 1 : NoB.total + NoB.input, 1) = u(NoB.generator + NoB.input + 1 : NoB.generator + 2*NoB.input, 1); %external phases

%% set dynamic state variables
if NoB.generator>0
    % generator phase angles
    delta = x(1:NoB.generator);
    % angular velocities
    omega = x(NoB.generator + 1 : 2*NoB.generator);
    % mechanical power
    P_m = x(2*NoB.generator + 1 : 3*NoB.generator);
    % commanded power production
    P_c = u(1:NoB.generator);
end

% set wind inputs
P_w(1 : NoB.total) = u(NoB.generator + 2*NoB.input + 1: NoB.generator + 2*NoB.input + NoB.total);

% algebraic variables
powVariables.E = E;
powVariables.V = V;
powVariables.Theta = Theta;
powVariables.P_w = P_w;

% dynamic variables (if generators exist)
if NoB.generator>0
    powVariables.delta = delta;
    powVariables.omega = omega;
    powVariables.P_m = P_m;
    powVariables.P_c = P_c;
end
    

% ------------------------------ END OF CODE ------------------------------
