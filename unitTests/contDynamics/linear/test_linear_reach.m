function res = test_linear_reach(~)
% test_linreach - unit_test_function of linear reachability analysis
%
% Checks the solution of the linearSys class for a small example
% with constant input against the analytical solution
%
% Syntax:  
%    res = test_linreach( ~ )
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
% Author:       Hendrik Roehm, Matthias Althoff
% Written:      02-March-2016
% Last update:  03-March-2016 (HR)
%               12-August-2016 (MA)
% Last revision:---

% Define acceptable reachability overapproximation error
% Since a linear system is tested this can be small
eps = 1e-10;

timeStep = 0.1;
numberOfTimeSteps = 10;

A = [0.1 1; -1 0.1]; % system matrix

%
% Configuration of the reachability analysis
%
options.x0 = [1; 1]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([0.5, 0.5])]); %initial set
options.timeStep = timeStep; %time step size in location 1
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.reductionTechnique = 'girard';
options.originContained = 0;
options.linAlg = 1;
%--------------------------------------------------------------------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %time step
    r = options.timeStep;
    %compute initial state factor
    options.factor(i)= r^(i)/factorial(i);    
end

%specify hybrid automaton--------------------------------------------------
%define large and small distance

B = zeros(2,1); % input matrix
linSys = linearSys('linearSys',A,B); %linear continuous system

%set input:
options.uTrans = 0; %input for simulation
options.U = zonotope(0); %input deviation from center

%simulate hybrid automaton
%linSys = simulate(linSys,options);

%compute reachable set
[Rnext, options] = initReach(linSys,options.R0,options);
for i=1:numberOfTimeSteps-1
    [Rnext, options] = post(linSys,[],options);
end

% Check solution against analytical solution by comparison of the vertices

% The solution can be computed analytically:
% x(t) = e^(A*t)*x0

% Vertices of both zonotopes
Vexact = get(vertices(expm(A)^(timeStep*numberOfTimeSteps)*options.R0), 'V');
Vcomputed = get(vertices(Rnext.tp), 'V');

% Sort both vertice lists so that corresponding vertices should be on the
% same index
Vexact = sortrows(Vexact')';
Vcomputed = sortrows(Vcomputed')';

% Each entry should be around zero;
res = all(all(abs(Vexact - Vcomputed) < eps));

%--------------------------------------------------------------------------


