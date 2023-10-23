function res = test_nonlinearSys_linError()
% test_nonlinearSys_linError - test if the linearization error for nonlinear
%                           systems is computed correctly
%
% Syntax:
%    res = test_nonlinearSys_linError()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       12-November-2018
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Parameters --------------------------------------------------------------

params.tFinal = 0.0021;
params.U = zonotope([0;0]);

x0 = zeros(18,1);
x0(4) = 19.9536;
x0(5) = 18.6195;
x0(8) = 0.7098;

G =   [[0 0 0 0.0203 0.0125 0.0123 0 0 0 0; ...
        0 0 0 -0.0068 -0.0015 0.0001 0.0071 0 0 0; ...
        0 0 0 0.1635 0.1524 -0.0041 0 0 0 0; ...
        0.15 0.1443 0.1391 0 0 0 0 0 0 0.06; ...
        -0.0016 -0.0045 -0.0072 0 0 0 0 0 0.0796 -0.0037; ...
        0 0 0 -0.0103 -0.0022 -0.0118 -0.0089 0.0786 0 0; ...
        0 0 0 0.0006 -0.0634 0 0 0 0 0; ...
        -7.1991 -0.0594 -0.0117 0 0 0 0 0 0 -0.0099]; eye(10)];

params.R0 = zonotope([x0,G]);


% Reachability Settings ---------------------------------------------------

options.timeStep = params.tFinal;     
options.taylorTerms = 5;
options.zonotopeOrder = 50;

options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@highorderBicycleDynamics);


% Reachability Analysis ---------------------------------------------------

% options check
options = validateOptions(sys,'reach',params,options);

% compute symbolic derivatives
derivatives(sys,options); 

% obtain factors for initial state and input solution time step
r = options.timeStep;
for i = 1:(options.taylorTerms+1)  
    options.factor(i) = r^(i)/factorial(i);    
end

% perform one reachability step
[R, options] = initReach(sys, options.R0, options);

% extract the set of linearization errors
err = R.tp{1}.error;
linError = zonotope([zeros(length(err),1),diag(err)]);

% evaluate the linearization error for a set of random points
p.u = center(params.U);
f0prev = sys.mFile(center(params.R0),p.u);
p.x = center(params.R0) + f0prev*0.5*options.timeStep;

f0 = sys.mFile(p.x,p.u);
[A,~] = sys.jacobian(p.x,p.u);

N = 5;
points = zeros(length(x0),N);

for i = 1:N
    p = randPoint(params.R0);

    temp = sys.mFile(p,[0;0]);
    points(:,i) = temp - A * (p - sys.linError.p.x) - f0;
end

% check if the set of linerization error contains all randomly computed points
linError = interval(linError);

linError = linError([1 3 5 6]);
points = points([1 3 5 6],:);

for i = 1:N
    if ~contains(linError,points(:,i))
    	res = false;
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
