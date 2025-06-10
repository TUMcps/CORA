function [OGain,tComp]= priv_observe_gain_PRadB(linsysDT,params,options)
% priv_observe_gain_PRadB - computes the gain for the guaranteed state estimation
%    approach from [1].
%
% Syntax:
%    [OGain,tComp] = priv_observe_gain_PRadB(linsysDT,params,options)
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
%    [1] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotope-based set-membership estimation for
%        multi-output uncertain systems. In Proc. of the IEEE
%        International Symposium on Intelligent Control (ISIC),
%        pages 212–217, 2013.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   02-January-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

timerVal = tic;

% E and F in [1] are chosen such that they are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
sys.E = generators(params.W);
sys.F = generators(params.V);

% obtain system dimension and nr of outputs
n = linsysDT.nrOfDims; 
nrOfOutputs = size(linsysDT.C,1); % not obj.nrOfOutputs since obj.C has been altered in the calling function!

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(n,n,'symmetric'); 
% gain matrix
Y = sdpvar(n,nrOfOutputs,'full'); 
% factor to be maximized for minimizing the P-radius
tau = sdpvar(1,1); 

% auxiliary value epsilon, above eq. (11) of [1]
phi = supremum(abs(interval(params.V)));
epsilon = norm(params.W) + phi'*phi;

%% optimization settings
beta_up = 1;  % Max value of beta
beta_lo = 0; % Min value of beta
beta_tol = 0.01;

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 1;  
options_sdp.shift = 1e-5;  

% Optimization loop
while(beta_up-beta_lo)> beta_tol
    beta_tst = (beta_up+beta_lo)/2;
    % implementation of symmetric matrix of eq. (17) in [1]
    SM = blkvar;
    SM(1,1) = beta_tst*P;
    SM(1,4) = linsysDT.A'*P-linsysDT.A'*linsysDT.C'*Y';
    SM(2,2) = sys.E'*sys.E;
    SM(2,4) = sys.E'*P-sys.E'*linsysDT.C'*Y';
    SM(3,3) = sys.F'*sys.F; 
    SM(3,4) = sys.F*Y'; % the paper uses Y'*sys.F, which is probably a mistake
    SM(4,4) = P;
    SM = sdpvar(SM);
    
    cond2 = (1-beta_tst)*P/epsilon;

    % define the optimization criterion as in eq. (17) of [1]
    % instead of maximizing tau, we minimize -tau
    objective = -tau; 

    % LMI problem to be solved
    constraint =  [(P>=0), (SM>=0), (cond2>=eye(n)*tau), (tau>=0)];

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
tComp = toc(timerVal);

% ------------------------------ END OF CODE ------------------------------
