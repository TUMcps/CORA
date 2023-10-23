function [OGain,tComp] = observe_gain_PRadE(obj,options)
% observe_gain_PRadE - computes the gain for the guaranteed state estimation
% approach from [1].
%
%
% Syntax:
%    [OGain,tComp] = observe_gain_PRadE(obj,options)
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
%    [1] Ye Wang, Meng Zhou, Vicenc Puig, Gabriela Cembrano, and
%        Zhenhua Wang. Zonotopic fault detection observer with H −
%        performance. In Proc. of the 36th IEEE Chinese Control
%        Conference, pages 7230–7235, 2017.
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

tic

% E and F in [1] are chosen such that they are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
sys.E = generators(options.W);
sys.F = generators(options.V);

% obtain system dimension and nr of outputs
n = obj.dim; 
nrOfOutputs = obj.nrOfOutputs;

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(n,n,'symmetric'); 
% gain matrix
Y = sdpvar(n,nrOfOutputs,'full'); 
        

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
    % implementation of symmetric matrix SM of eq. (17) in [1]
    SM = blkvar;
    SM(1,1) = -beta_tst*P;
    SM(1,4) = obj.A'*P-obj.C'*Y';
    SM(2,2) = -sys.E'*sys.E;
    SM(2,4)= sys.E'*P;
    SM(3,3) = -sys.F'*sys.F; 
    SM(3,4) = sys.F*Y';
    SM(4,4) = -P;
    SM = sdpvar(SM);

    % optimization criterion seems to not be specified in [1]
    objective = -trace(P); 

    % LMI problem to be solved
    constraint = [(P>=0), (SM<=0)];

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
