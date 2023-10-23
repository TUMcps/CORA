function [OGain,tComp]= observe_gain_PRadD(obj,options)
% observe_gain_PRadD - computes the gain for the guaranteed state estimation
% approach from [1].
%
% Syntax:
%    [OGain,tComp]= observe_gain_PRadD(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    OGain - observer gain
%    tComp - computation time
%
% Reference:
%    [1] Ye Wang, Zhenhua Wang, Vicenc Puig, and Gabriela
%        Cembrano. Zonotopic set-membership state estimation for
%        discrete-time descriptor LPV systems. IEEE Transactions
%        on Automatic Control, 64(5):2092-2099, 2019.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   01-March-2021
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

% choice of alpha and beta depends on designer
alpha = 0.5; % has to be in ]0,1[
beta = 0.5; % has to be in ]0,1[

% define required matrices
E_new = [E, zeros(size(E,1), size(F,2))];
F_new = [zeros(size(F,1), size(E,2)), F];

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(n,n,'symmetric'); 
% gain matrix
Y = sdpvar(n,nrOfOutputs,'full'); 
% identity matrix
I = eye(size(E_new,2));
I2 = eye(size(F_new,2));
        

%% optimization settings
gamma_up = 100;  % Max value of beta
gamma_lo = 0; % Min value of beta
gamma_tol = 0.01;

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver; 
options_sdp.shift = 1e-5;
options_sdp.verbose = 1;  

% Optimization loop
while(gamma_up-gamma_lo)> gamma_tol
    gamma_tst = (gamma_up+gamma_lo)/2;
    % implementation of the symmetric matrix SM in eq. (33) in [1]
    % choice of T and N in (4) of [1]: T = I, N = 0
    SM = blkvar;
    SM(1,1) = alpha*P;
    SM(2,2) = (1-alpha)*beta*I;
    SM(3,3) = (1-alpha)*(1-beta)*I2;
    SM(4,1) = (P-Y*obj.C)*obj.A;
    SM(4,2) = (P-Y*obj.C)*E_new;
    SM(4,3) = Y*F_new;
    SM(4,4) = P;
    SM = sdpvar(SM);
    
    % implementation of the symmetric matrix SM in eq. (20) in [1]
    SM2 = blkvar;
    SM2(1,1) = eye(n);
    SM2(2,1) = gamma_tst*P;
    SM2(2,2) = P;
    SM2 = sdpvar(SM2);

    % optimization criterion seems to not be specified in [1]
    objective = -trace(P); 

    % LMI problem to be solved
    constraint =  [(P>=0), (SM>=0), (SM2<=0)];

    % Solve LMI conditions
    solpb = optimize(constraint, objective, options_sdp);

    % Check if LMI is feasible
    if solpb.problem == 1
        disp('LMIs are infeasible');
        gamma_up = gamma_tst;
    else
        gamma_lo = gamma_tst;
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
