function [OGain,tComp]= observe_gain_NomG(obj,options)
% observe_gain_NomG - computes the gain for the guaranted state estimation
% approach from [1].
%
% Syntax:  
%    [R,Rout] = observe_NomG(obj,options)
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
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       05-Jan-2021
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------

tic

% obtain system dimension and nr of outputs
dim = size(obj.A,1); 
nrOfOutputs = size(obj.C,1);

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 0;  

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(dim,dim,'symmetric'); 
% gain matrix
Y = sdpvar(dim,nrOfOutputs,'full'); 
        
mu =1; % decay rate

%% Solving the LMI Problem
% create symmetric matrix SM of LMI
SM = blkvar;
SM(1,1) = mu*P;
SM(2,1) = P*obj.A-Y*obj.C;
SM(2,2) = P;
SM = sdpvar(SM);
ObjR = [P>=0,SM>=0]; % The objective function
sol_lmi = optimize(ObjR,[],options_sdp); % optimze the LMIs

% Check if LMI is feasible
if sol_lmi.problem == 1
    disp('LMIs are infeasible');
    OGain= zeros(dim,nrOfOutputs);
else
    % extract values from YALMIP symbolic decision variables
    P = value(P);
    Y = value(Y);

    % store gain
    OGain = P\Y;
end

% computation time
tComp = toc;

%------------- END OF CODE --------------