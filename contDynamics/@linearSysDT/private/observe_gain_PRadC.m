function [OGain,tComp]= observe_gain_PRadC(obj,options)
% observe_gain_PRadC - computes the gain for the guaranteed state estimation
%    approach from [1] and [2].
%
% Syntax:
%    [OGain,tComp] = observe_gain_PRadC(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    OGain - observer gain
%    tComp - computation time
%
% Example: 
%    -
%
% Reference:
%    [1] Ye Wang, Teodoro Alamo, Vicenc Puig, and Gabriela
%        Cembrano. A distributed set-membership approach based on
%        zonotopes for interconnected systems. In Proc. of the IEEE
%        Conference on Decision and Control (CDC), pages 668–673, 2018.
%    [2] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
%    [3] Ye Wang, Zhenhua Wang, Vicenc Puig, and Gabriela
%        Cembrano. Zonotopic set-membership state estimation for
%        discrete-time descriptor LPV systems. IEEE Transactions
%        on Automatic Control, 64(5):2092-2099, 2019.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   02-January-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic;

% E and F in [1] are chosen such that they are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
E = generators(options.W);
F = generators(options.V);

% obtain system dimension and nr of outputs
n = obj.dim; 
nrOfOutputs = obj.nrOfOutputs;

%% define YALMIPs symbolic decision variables
% weight of F_P radius
P = sdpvar(n,n,'symmetric'); 
% gain matrix
Y = sdpvar(n,nrOfOutputs,'full'); 

%% set G and O 
G = sdpvar(size(E,2),size(E,2),'diag');
O = sdpvar(nrOfOutputs,nrOfOutputs,'diag');     

%% optimization settings
beta_up = 1;  % Max value of beta
beta_lo = 0; % Min value of beta
beta_tol = 0.01;

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver; 
options_sdp.shift = 1e-5;
options_sdp.verbose = 1;  

% Optimization loop
while(beta_up-beta_lo)> beta_tol
    beta_tst = (beta_up+beta_lo)/2;
    % choice of T and N in (4) of [3]: T = I, N = 0
    % implementation of symmetric matrix SM of eq. (18) in [2]
    % the 3rd row and column can be removed due to N = 0
    % W --> P
    % \gamma --> \beta
    % \Gamma --> G
    % E --> I (identity; only required for descriptor systems)
    %
    % this is also presented in [1] on p. 671
    
    SM = blkvar;
    SM(1,1) = beta_tst*P;
    SM(2,2) = G;
    SM(3,3) = O;
    SM(4,1) = (P-Y*obj.C)*obj.A;
    SM(4,2) = (P-Y*obj.C)*E;
    SM(4,3) = Y*F;
    SM(4,4) = P;
    SM = sdpvar(SM);

    % optimization criterion seems to not be specified in [1]
    objective = -trace(P); 

    % LMI problem to be solved
    constraint =  [(P>=0), (SM>=0) ,(G>=0), (O>=0)];

    % Solve LMI conditions
    solpb = optimize(constraint, objective, options_sdp);

    % Check if LMI is feasible
    if solpb.problem == 1
        disp('LMIs are infeasible');
        beta_up = beta_tst;
    else
        beta_lo = beta_tst;
    end
end
% extract values from YALMIP symbolic decision variables
P = value(P);
Y = value(Y);

% store gain
OGain = P\Y;

% computation time
tComp = toc;

% ------------------------------ END OF CODE ------------------------------
