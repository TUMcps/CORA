function generateConstraintEquations(scenario,powVariables)
% generateConstraintEquations - generates the constraint equations of a 
%    power system, more information can be found in [1, Sec. VII].
%
% Syntax:
%    generateConstraintEquations(scenario,powVariables)
%
% Inputs:
%    scenario - struct specifying a power system scenario
%    powVariables - symbolic variables for the power system
%
% Outputs:
%    -
%
% References:
%    [1] M. Althoff, "Benchmarks for the Formal Verification of Power 
%        Systems", Proc. of the 9th International Workshop on Applied 
%        Verification of Continuous and Hybrid Systems, 
%        2022, x-x

% Authors:       Matthias Althoff
% Written:       15-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain certain structs and variables
bus = scenario.bus;
genParam = scenario.genParam;
Ycomplex = scenario.Y;
Pd = scenario.Pd;
Qd = scenario.Qd;

% number of buses (NoB)
NoB.input = length(bus.input);
NoB.generator = length(bus.generator);
NoB.load = length(bus.load);
NoB.total = length(scenario.Y) - NoB.input;

% remove connection to generator
if ~isempty(bus.fault)
    genParam.Y_m(bus.fault) = 0;
end

% obatain absolute values and angles
Y = abs(Ycomplex);
delta_load = angle(Ycomplex);

% obatain power system variables
V = powVariables.V;
Theta = powVariables.Theta;
P_w = powVariables.P_w;
if NoB.generator>0
    E = powVariables.E;
    delta = powVariables.delta;
end

% set P and Q for different buses
% generator buses; see equation above eq. (2) in [1]
for i = 1:NoB.generator
    if genParam.Psi_g(i) == -pi/2 % use special case to simplify later symbolic derivations
        P_g(i) = genParam.Y_m(i)*E(i)*V(i)*sin(delta(i) - Theta(i));
        Q_g(i) = genParam.Y_m(i)*E(i)*V(i)*cos(delta(i) - Theta(i)) - genParam.Y_m(i)*V(i)^2 ;
    else
        P_g(i) = genParam.Y_m(i)*E(i)*V(i)*cos(genParam.Psi_g(i) + delta(i) - Theta(i)) - genParam.Y_m(i)*V(i)^2*cos(genParam.Psi_g(i));
        Q_g(i) = - genParam.Y_m(i)*E(i)*V(i)*sin(genParam.Psi_g(i) + delta(i) - Theta(i)) + genParam.Y_m(i)*V(i)^2*sin(genParam.Psi_g(i));
    end
    % combine powers
    P(i) = P_g(i) - Pd(i) + P_w(i);
    Q(i) = Q_g(i) - Qd(i);
end
% load buses
for i = (NoB.generator+1) : NoB.total
    P(i) = -Pd(i) + P_w(i); 
    Q(i) = -Qd(i);
end


% the first set of symbolic variables represents active power, 
% the second set reactive power

% active power; eq. (3) in [1]
for i = 1:NoB.total
    %init
    g(i, 1) = sym(0);
    for n = 1:(NoB.total + NoB.input)
        g(i, 1) = g(i, 1) + V(i)*V(n)*Y(i,n)*cos(delta_load(i,n) + Theta(n) - Theta(i));
    end
    g(i, 1) = g(i, 1) - P(i);
end
    
% reactive power; eq. (3) in [1]
for i = 1:NoB.total
    %init
    g(NoB.total + i, 1) = sym(0);
    for n = 1:(NoB.total + NoB.input)
        g(NoB.total + i, 1) = g(NoB.total + i, 1) - (V(i)*V(n)*Y(i,n)*sin(delta_load(i,n) + Theta(n) - Theta(i)));
    end
    g(NoB.total + i, 1) = g(NoB.total + i, 1) - Q(i);
end

% consider slack bus
if ~isempty(bus.slack)
    %add first active and reactive power equation together
    g(1, 1) = g(NoB.total + 1, 1) + g(1, 1);
    %remove first active power equation
    g(NoB.total + 1) = [];
end

% create file
createFileFromFunction(g,[scenario.name,'_con'],'g','x,y,u');

% ------------------------------ END OF CODE ------------------------------
