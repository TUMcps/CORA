function scenarioNew = reorderBuses(scenario)
% reorderBuses - rename buses such that 
% 1 is the slack bus
% 2...(nrOfGenBuses+1) are the generator buses
% (nrOfGenBuses+2)...(nrOfGenBuses+nrOfLoadBuses+1) are the load buses
% (nrOfGenBuses+nrOfLoadBuses+2)...(nrOfGenBuses+nrOfLoadBuses+nrOfInputBuses+1) are the input buses
%
%
% More information can be found in [1, Sec. VII].
%
% Syntax:
%    scenarioNew = reorderBuses(scenario)
%
% Inputs:
%    scenario - struct specifying a power system scenario
%
% Outputs:
%    scenarioNew - struct specifying a power system scenario
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

% copy old values
scenarioNew = scenario;

% remove entries so that dimensions match the ones of the subsystems
scenarioNew.Y = [];
scenarioNew.Pd = [];
scenarioNew.Qd = [];
scenarioNew.VM = [];

% create order vector
bus = scenario.bus;
orderVec = [bus.generator'; bus.load'; bus.input'];

%% reorder admittance matrix
for i=1:length(orderVec)
    for j=1:length(orderVec)
        %get new indices
        i_old = orderVec(i);
        j_old = orderVec(j);
        
        %write value to new indices
        scenarioNew.Y(i, j) = scenario.Y(i_old,j_old);
    end
end

%% reorder demanded power
for i=1:length(orderVec)
    %get new indices
    i_old = orderVec(i);
    
    %write value to new indices
    scenarioNew.Pd(i) = scenario.Pd(i_old);
    scenarioNew.Qd(i) = scenario.Qd(i_old);
end

%% reorder bus voltages
for i=1:length(bus.generator)
    %get new indices
    i_old = bus.generator(i);
    
    %write value to new indices
    scenarioNew.VM(i) = scenario.VM(i_old);
end

% ------------------------------ END OF CODE ------------------------------
