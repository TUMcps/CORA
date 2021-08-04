function [OGain,tComp]= observe_gain_HinfG(obj,options)
% observe_gain_HinfG - computes the gain for the guaranted state estimation
% approach from [1].
%
% Syntax:  
%    [OGain,tComp]= observe_gain_HinfG(obj,options)
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
%    [1] W. Tang, Z. Wang, Y. Wang, T. Raissi, and Y. Shen.
%        Interval estimation methods for discrete-time linear time-
%        invariant systems. IEEE Transactions on Automatic Control,
%        64(11):4717-4724, 2019.
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

% It is assumed that E and F are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
E = generators(options.W);
F = generators(options.V);

% obtain system dimension and nr of outputs
dim = size(obj.A,1); 
nrOfOutputs = size(obj.C,1);
nw = size(E,2);
nv = size(F,1);

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 0;  

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(dim,dim,'symmetric'); 
% gain matrix
Y = sdpvar(dim,nrOfOutputs,'full'); 
        
gamma = linspace(1,50,100); % todo: description

%% Solving the LMI Problem
for i = 1:length(gamma)
    % compute symmetric matrix SM of LMI
    SM= blkvar;
    SM(1,1) = eye(dim)-P;
    SM(2,1) = 0;
    SM(2,2) = -gamma(i)*gamma(i)*eye(nw);
    SM(3,1) = 0;
    SM(3,2) = 0;
    SM(3,3) = -gamma(i)*gamma(i)*eye(nv);
    SM(4,1) = P*obj.A-Y*obj.C;
    SM(4,2) = P*E;
    SM(4,3) = -Y*F;
    SM(4,4) = -P;
    SM= sdpvar(SM);
    ObjR= [P>=0,SM<=0,gamma(i)>=0]; % The objective function    
    solpb = optimize(ObjR,gamma(i)*gamma(i),options_sdp); % optimze the LMIs
    
    % Check if LMI is feasible
    if solpb.problem == 1
        disp('LMIs are infeasible');
        OGain= zeros(dim,nrOfOutputs);
    else
        % extract values from YALMIP symbolic decision variables
        P = value(P);
        Y = value(Y);

        % store gain
        OGain = P\Y;
        break
    end
end

% computation time
tComp = toc;

 
end

%------------- END OF CODE --------------