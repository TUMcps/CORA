function index = indexForSubsystems(sys, subsystem, bus, bus_subsystem, statesPerBus)
% indexForSubsystems - relates states and inputs of a system to its
% subsystems; currently this is only implemented for power systems.
%
% The following symbolic variables are used: 
% 1 is the slack bus
% 2...(nrOfGenBuses+1) are the generator buses
% (nrOfGenBuses+2)...(nrOfGenBuses+nrOfLoadBuses+2) are the load buses
%
% dynamic states:
% x(i - 1) = delta(i), i = 2...nrOfGenerators; generator phase angles
% x(nrOfGenerators - 1 + i) = omega(i), i = 1...nrOfGenerators; angular velocity
% x(2*nrOfGenerators - 1 + i) = T_m(i), i = 1...nrOfGenerators; torque
%
% algebraic states:
% y(i) = E(i), i = 1...nrOfGenerators; generator voltage
% y(nrOfGenerators + i) = V(nrOfGenerators + i), i = 1...nrOfLoadBuses; bus voltage
% y(nrOfGenerators + nrOfLoadBuses + i) = Theta(i) - delta, j = 1...nrOfBuses
%
% inputs:
% u(i) = P_c(i), i = 1...nrOfGenerators; commanded power production
% u(nrOfGenerators + i) = P^d_g(i), i = 1...nrOfBuses; directly injected active power
% u(nrOfGenerators + nrOfBuses + i) = Q^d_g(i), i = 1...nrOfBuses; directly injected reactive power
% u(nrOfGenerators + 2*nrOfBuses + i) = V(i), i = 1...nrOfInputs; input voltages
% u(nrOfGenerators + 2*nrOfBuses + nrOfInputs + i) = Theta(i), i = 1...nrOfInputs; input phases
%
%
% Syntax:
%    index = indexForSubsystems(sys, subsystem, bus, bus_subsystem, statesPerBus)
%
% Inputs:
%    sys - full system
%    subsystem - cell array of subsystems
%    bus - information of buses of the full system
%    bus_subsystem - information of buses of the subsystems
%    statesPerBus - nr of continuous state variables per bus
%
% Outputs:
%    index - index struct relating states and inputs of a system to its
%            subsystems
%
% Reference:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% returns the inputs to the goalSubsystem; the function also works for sets
% as inputs

% global quantities
globalNrOfGenerators = length(bus.generator);
globalNrOfBuses = length(bus.generator) + length(bus.load);
globalNrOfInputBuses = length(bus.input);

% create all buses
bus.all = [bus.generator'; bus.load'];

% local quantities
for iSys = 1:length(subsystem)
    nrOfGenerators{iSys} = length(bus_subsystem{iSys}.generator);
    nrOfBuses{iSys} = length(bus_subsystem{iSys}.generator) + length(bus_subsystem{iSys}.load);
    nrOfInputBuses{iSys} = length(bus_subsystem{iSys}.input);
    
    % create all buses
    bus_subsystem{iSys}.all = [bus_subsystem{iSys}.generator'; bus_subsystem{iSys}.load'];
end


% check all subsystems
for iSys = 1:length(subsystem)
    % init index matrices
    indX = zeros(subsystem{iSys}.dim, sys.dim); %index for x-values
    indY = zeros(subsystem{iSys}.nrOfConstraints, sys.nrOfConstraints); %index for y-values
    indUy = zeros(subsystem{iSys}.nrOfInputs, sys.nrOfConstraints); %index for u-values
    indUv = zeros(subsystem{iSys}.nrOfInputs, globalNrOfGenerators); %index for u-values from constant voltage
    indUu = zeros(subsystem{iSys}.nrOfInputs, sys.nrOfInputs); %index for u-values from inputs
      
    %% obtain dynamic state mapping and parts of indUu
    % check all buses within the subsystem
    for i = 1:nrOfGenerators{iSys}
        % obtain bus of current subsystem
        busNr = bus_subsystem{iSys}.generator(i);
        % obtain global bus index
        globalBus = find(bus.all == busNr);
        % build indX
        for iState = 1:statesPerBus
            indState = globalNrOfGenerators*(iState-1) + globalBus;
            indStateNew = nrOfGenerators{iSys}*(iState-1) + i;
            indX(indStateNew, indState) = 1;
        end
        % build indUu; the first entries are the commanded power production
        indState = globalBus;
        indStateNew = i;
        indUu(indStateNew, indState) = 1;
    end
    
    %% obtain algebraic state mapping and parts of indUu
    % check all buses within the subsystem
    for i = 1:nrOfBuses{iSys}
        % obtain bus of current subsystem
        busNr = bus_subsystem{iSys}.all(i);
        % obtain global bus index
        globalBus = find(bus.all == busNr);
        % build indY
        if globalBus ~= 1 % global bus is not the slack bus
            indState1 = globalBus; 
            indState2 = globalNrOfBuses + globalBus - 1;
            indStateNew1 = i;
            if bus_subsystem{iSys}.slack
                indStateNew2 = nrOfBuses{iSys} + i - 1;
            else
                indStateNew2 = nrOfBuses{iSys} + i;
            end
            indY(indStateNew1, indState1) = 1;
            indY(indStateNew2, indState2) = 1;
        else % global bus is the slack bus
            indY(1, 1) = 1;
        end
        % build indUu; requires further checking
        indState = globalNrOfGenerators + 2*globalNrOfInputBuses + globalBus;
        indStateNew = nrOfGenerators{iSys} + 2*nrOfInputBuses{iSys} + i;
        indUu(indStateNew, indState) = 1;
    end
    
    %% obtain mapping from algebraic states to inputs
    % check all buses within the subsystem
    for i = 1:nrOfInputBuses{iSys}
        %obtain bus of current subsystem
        busNr = bus_subsystem{iSys}.input(i);
        %obtain global bus index
        globalBus = find(bus.all == busNr);
        %check if generator bus
        isGeneratorBus = ~isempty(find(bus.generator == busNr, 1));
        
        %build indUy, indUv
        if globalBus ~= 1 %global bus is not the slack bus

            indState1 = globalBus; 
            indState2 = globalNrOfBuses + globalBus - 1;
            indStateNew1 = nrOfGenerators{iSys} + i;
            indStateNew2 = nrOfGenerators{iSys} + nrOfInputBuses{iSys} + i;
            
            %generator bus
            if isGeneratorBus
                indUv(indStateNew1, indState1) = 1;
                indUy(indStateNew2, indState2) = 1;
            else
                indUy(indStateNew1, indState1) = 1;
                indUy(indStateNew2, indState2) = 1;
            end
            
        else
            indUv( nrOfGenerators{iSys} + 1, 1) = 1;
        end
    end      
    
    %store values
    index{iSys}.X = indX;
    index{iSys}.Y = indY;
    index{iSys}.Uv = indUv;
    index{iSys}.Uy = indUy;
    index{iSys}.Uu = indUu;
end

% ------------------------------ END OF CODE ------------------------------
