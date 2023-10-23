function generateDynamicEquations(scenario, powVariables)
% generateDynamicEquations - generates the dynamic equations of a power 
%    system, more information can be found in [1, Sec. VII].
%
% Syntax:
%    generateDynamicEquations(scenario, powVariables)
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

% obtain bus and genParam struct
bus = scenario.bus;
genParam = scenario.genParam;

% number of buses (NoB)
NoB.generator = length(bus.generator);

if NoB.generator>0
    %remove connection to generator
    if ~isempty(bus.fault)
        genParam.Y_m(bus.fault) = 0;
    end

    %obatain power system variables
    E = powVariables.E;
    V = powVariables.V;
    Theta = powVariables.Theta;
    delta = powVariables.delta;
    omega = powVariables.omega;
    P_m = powVariables.P_m;
    P_c = powVariables.P_c;


    %the first set of symbolic variables represents relative angles, the second 
    %set angular velocities, the third set torques

    %relative angles; first angle is reference and not specified; eq. (18) in [1]
    for i = 1:NoB.generator
        f(i, 1) = omega(i) - genParam.omega_s;
    end

    % angular velocities; eq. (1) in [1]
    for i = 1:NoB.generator
        % obtain P_g
        if genParam.Psi_g(i) == -pi/2 % use special case to simplify later symbolic derivations
            P_g(i) = genParam.Y_m(i)*E(i)*V(i)*sin(delta(i) - Theta(i));
        else
            P_g(i) = genParam.Y_m(i)*E(i)*V(i)*cos(genParam.Psi_g(i) + delta(i) - Theta(i)) - genParam.Y_m(i)*V(i)^2*cos(genParam.Psi_g(i));
        end
        % right-hand side of dynamic equations
        f(NoB.generator + i, 1) = -genParam.D(i)/genParam.M(i)*(omega(i) - genParam.omega_s)...
            + 1/genParam.M(i)*P_m(i) - 1/genParam.M(i)*P_g(i);
    end

    % mechanical power; eq. (1) in [1]
    for i = 1:NoB.generator
        %init
        f(2*NoB.generator + i, 1) = -1/(genParam.T_sv(i)*genParam.R_d(i))*(omega(i) - genParam.omega_s)...
            - 1/genParam.T_sv(i)*P_m(i) + 1/genParam.T_sv(i)*P_c(i);
    end
else
    f = [];
end

% create file
createFileFromFunction(f,[scenario.name,'_dyn'],'f','x,y,u');


% ------------------------------ END OF CODE ------------------------------
