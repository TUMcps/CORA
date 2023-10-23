function [R,tcomp] = observe_ESO_D(obj,options)
% observe_ESO_D - computes the guaranteed state estimation approach of [1].
%
%
% Syntax:
%    [R,tcomp] = observe_ESO_D(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] John J. Martinez, Nassim Loukkas, and Nacim Meslem.
%        H∞ set-membership observer design for discrete-time LPV
%        systems. International Journal of Control, 93(10):2314-2325,
%        2020.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains

[L,P,gamma,lambda] = observe_gain_ESO_D(obj,options);

%%initialize computation
%time period
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;
timeSteps = length(tVec);

% initialize parameter for the output equation
R = cell(length(tVec),1);

% store first reachable set
Rnext.tp = options.R0;
R{1} = Rnext.tp;
x = Rnext.tp.center;

% compute e_inf, eq. (69) in [1]
n = obj.dim; 
nrOfOutputs = size(obj.C,1);
w_bar = ones(n+nrOfOutputs,1); % without loss of generality, the disturbances are bounded by 1
c_bar = gamma^2/lambda*(w_bar'*w_bar);
e_inf = diag((P/c_bar)^-0.5);

% init mu according to (74) in [1]
% enclosing ball of initial set
r = radius(options.R0);
lambda_max = eigs(P,1);
mu = lambda_max*r^2/c_bar;

tic;

%% loop over all time steps
for k = 1:timeSteps-1
    
    % center, eq. (65) in [1]
    x = obj.A*x + obj.B*options.uTransVec(:,k) + obj.c + L*(options.y(:,k) - obj.C*x);
    % update rho, eq. (66) in [1]
    mu = (1-lambda)*mu + lambda;
    
    % upper and lower bound, eq. (67), (68) in [1]
    delta = e_inf*mu^0.5;
    Rnext.tp = interval(x - delta, x + delta);

    % Store result
    R{k+1} = Rnext.tp;
end


tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
