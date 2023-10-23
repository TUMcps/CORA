function res = test_linearSys_reach_01
% test_linearSys_reach_01 - unit test function of linear reachability analysis
%
% Checks the solution of the linearSys class for a small example
% with constant input against the analytical solution
%
% Syntax:
%    res = test_linearSys_reach_01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Hendrik Roehm, Matthias Althoff
% Written:       02-March-2016
% Last update:   03-March-2016 (HR)
%                12-August-2016 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Define acceptable reachability overapproximation error
% Since a linear system is tested this can be small
eps = 1e-10;

timeStep = 0.1;
numberOfTimeSteps = 50;

A = [0.1 1; -1 0.1]; % system matrix
B = 1; % input matrix
linSys = linearSys('linearSys',A,B); %linear continuous system

% Configuration of the reachability analysis ------------------------------
options.R0 = zonotope([[1; 1], diag([0.5, 0.5])]); %initial set
options.timeStep = timeStep;
options.tFinal = options.timeStep * numberOfTimeSteps;
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.compTimePoint = true;
% note: warning for compTimePoint is ok, since post_Euclidean used
options.U = zonotope(0);
options.uTrans = 0;
options.originContained = true;
options.reductionTechnique = 'girard';
options.linAlg = 'standard';
%--------------------------------------------------------------------------

%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

%compute reachable set
[Rnext, options] = initReach(linSys,options.R0,options);
for i=1:numberOfTimeSteps-1
    [Rnext, options] = post(linSys,[],options);
end

% Check solution against analytical solution by comparison of the vertices

% The solution can be computed analytically:
% x(t) = e^(A*t)*x0

% Vertices of both zonotopes
Vexact = vertices(expm(A)^(timeStep*numberOfTimeSteps)*options.R0);
Vcomputed = vertices(Rnext.tp);

% Sort both vertice lists so that corresponding vertices should be on the
% same index
Vexact = sortrows(Vexact')';
Vcomputed = sortrows(Vcomputed')';

% Each entry should be around zero;
res = all(all(abs(Vexact - Vcomputed) < eps));


% ------------------------------ END OF CODE ------------------------------
