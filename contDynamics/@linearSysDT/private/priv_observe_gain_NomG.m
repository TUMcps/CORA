function [OGain,tComp]= priv_observe_gain_NomG(linsysDT,params,options)
% priv_observe_gain_NomG - computes the gain for the guaranteed state estimation
%    approach from [1].
%
% Syntax:
%    [R,Rout] = priv_observe_gain_NomG(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    OGain - observer gain
%    tComp - computation time
%
% Reference:
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic;

% obtain system dimension and nr of outputs
n = linsysDT.nrOfDims;
nrOfOutputs = linsysDT.nrOfOutputs;

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 0;  

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(n,n,'symmetric'); 
% gain matrix
Y = sdpvar(n,nrOfOutputs,'full'); 
        
mu = 1; % decay rate

%% Solving the LMI Problem
% create symmetric matrix SM of LMI
SM = blkvar;
SM(1,1) = mu*P;
SM(2,1) = P*linsysDT.A-Y*linsysDT.C;
SM(2,2) = P;
SM = sdpvar(SM);
constraint = [P>=0,SM>=0]; % constraint
sol_lmi = optimize(constraint, [], options_sdp); % solve LMIs; no objective

% Check if LMI is feasible
if sol_lmi.problem == 1
    disp('LMIs are infeasible');
else
    % extract values from YALMIP symbolic decision variables
    P = value(P);
    Y = value(Y);

    % store gain
    OGain = P\Y;
end

% computation time
tComp = toc;

% ------------------------------ END OF CODE ------------------------------
