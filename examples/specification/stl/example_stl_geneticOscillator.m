function res = example_stl_geneticOscillator
% example_stl_geneticOscillator - example of signal temporal
% logic checking of the genetic oscillator model from [1]
%
% Syntax:
%    res = example_stl_geneticOscillator()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean
%
% References:
%    [1] T. Wright, I. Stark. "Property-directed verified monitoring of
%        signal temporal logic", Runtime Verification 2020

% Authors:       Benedikt Seidl
% Written:       04-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Warning
disp('Warning: This example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

% Parameter ---------------------------------------------------------------

params.R0 = zonotope(interval( ...
    [0.98; 1.28; 0.08; 0.08; 0.08; 1.28; 2.48; 0.58; 1.28], ...
    [1.02; 1.32; 0.12; 0.12; 0.12; 1.32; 2.52; 0.62; 1.32]));
params.tFinal = 4.5;


% Reachability Settings ---------------------------------------------------

options.alg = 'poly-adaptive';


% System Dynamics ---------------------------------------------------------

f = @(x,u) [50 * x(3) - 0.1 * x(1) * x(6);
    100 * x(4) - x(1) * x(2);
    0.1 * x(1) * x(6) - 50 * x(3);
    x(2) * x(6) - 100 * x(4);
    5 * x(3) + 0.5 * x(1) - 10 * x(5);
    50 * x(5) + 50 * x(3) + 100 * x(4) - x(6) * (0.1 * x(1) + x(2) + 2 * x(8) + 1);
    50 * x(4) + 0.01 * x(2) - 0.5 * x(7);
    0.5 * x(7) - 2 * x(6) * x(8) + x(9) - 0.2 * x(8);
    2 * x(6) * x(8) - x(9)];

sys = nonlinearSys(f);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

tic
simOpt.points = 5;
simRes = simulateRandom(sys, params, simOpt);
tComp = toc;
disp(['computation time of simulated traces: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

% 0.032 − 125^2(x4 − 0.003)^2 − 3(x6 − 0.5)^2 > 0
% 125^2(x4 - 0.003)^2 + 3(x6 - 0.5)^2 < 0.032
% 488281.25(x4 - 0.003)^2 + 93.75(x6 - 0.5)^2 < 1

q = [0.003; 0.5];
Q = [1/488281.25 0;
        0 1/93.75];

factor = 1.5;

% define atomic propositions

% x6 - 1 > 0
P = stl('P', atomicProposition(halfspace([0 0 0 0 0 -1 0 0 0], -1)));

% 0.032 − 125^2(x4 − 0.003)^2 − 3(x6 − 0.5)^2 > 0
Q = stl('Q', atomicProposition(ellipsoid(Q * factor, q), [4 6]));

% define formula
phi = globally(P | globally(Q, interval(3,3.5)), interval(0,1));


res = true;

% Model check formula on reachable set
tic
res = res && modelChecking(R,phi,'signals');
tComp = toc;
disp(['verification of reachable set: ',num2str(tComp)]);


% Verify all simulation traces
tic
res = res && monitorSTL(simRes,phi);
tComp = toc;
disp(['verification of simulation traces: ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

dims = {[4 6]};

for k = 1:length(dims)
    % start figure
    figure; hold on; box on;
    projDim = dims{k};

    % plot reachable sets
    plot(R,projDim);

    % plot initial set
    plot(params.R0,projDim,'k','FaceColor','w');

    % plot simulation results
    plot(simRes,projDim);

    % plot atomic propositions
    plot(polytope(P.lhs.set), projDim, 'y');
    plot(Q.lhs.set, [1 2], 'g');

    % label plot
    xlabel(['x_{',num2str(projDim(1)),'}']);
    ylabel(['x_{',num2str(projDim(2)),'}']);
end

end

% ------------------------------ END OF CODE ------------------------------
