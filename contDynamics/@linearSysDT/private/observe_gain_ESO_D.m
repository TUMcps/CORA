function [OGain,P,gamma,lambda,tComp] = observe_gain_ESO_D(obj,options)
% observe_gain_ESO_D - computes the gain for the guaranted state estimation
% approach from [1]. 
%
%
% Syntax:  
%    [OGain,P,gamma,lambda,tComp] = observe_gain_ESO_D(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    OGain - observer gain
%    P - result of Alg. 1 in [1]
%    gamma - result of Alg. 1 in [1]
%    lambda - minimum generalized eigenvalue, see line 7 of Alg. 1 in [1]
%    tComp - computation time
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

% Author:        Matthias Althoff
% Written:       04-Mar-2021
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------

tic

% essential values
dim = size(obj.A,1); 
nrOfOutputs = size(obj.C,1);

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 0;  

%% Alg. 1 of [1] 
% initialization of Alg. 1
Q = eye(length(obj.A)); 
% solve LMI, line 1 of Alg. 1
[P,Y] = solveLMI(obj,Q,options,options_sdp);
% compute L, line 2 of Alg. 1
L = P\Y;
% Compute observer matrix, line 3 of Alg. 1
Ao = obj.A-L*obj.C;
% Compute E, line 4 of Alg. 1
E = [options.W.Q, -L*options.V.Q];
% Estimate the expected steady-state covariance matrix V, line 5 of Alg. 1
D_w = diag(ones(1,dim+nrOfOutputs)); % without loss of generality, the disturbances are bounded by 1
W = (1/3)*D_w; % eq. (64)
V = sdpvar(dim,dim,'symmetric');
Fa = [V>=0, Ao*V*Ao'-V<=-E*W*E']; % eq. (63)
optimize(Fa);    
V = value(V);
% Obtain matrix Q, line 6 of Alg. 1
Q = inv(V);
% solve LMI, line 7 of Alg. 1
[P,Y,gamma] = solveLMI(obj,Q,options,options_sdp);
% compute gain, line 8 of Alg. 1
OGain = P\Y;
% Compute lambda as the minimum generalised eigenvalue of the pair (Q,P), line 9 of Alg. 1
lambda = min(eigs(Q,P));


% computation time
tComp = toc;
end


function [P,Y,gamma] = solveLMI(obj,Q,options,options_sdp)
    
% shape matrices of the disturbance and porcess noise sets
F = options.W.Q;
E = options.V.Q;

% obtain system dimension and nr of outputs
dim = size(obj.A,1); 
nrOfOutputs = size(obj.C,1);
nrOfDistGens = size(F,2);

%% define YALMIPs symbolic decision variables
% state
P = sdpvar(dim,dim,'symmetric'); 
% gain matrix
Y = sdpvar(dim,nrOfOutputs,'full'); 
% identity matrix
I = eye(nrOfDistGens + nrOfOutputs);

% possible values of gam
gam = linspace(1,50,100);

% Optimization loop
for k = 1:length(gam)
    
    % implementation of symmetric matrix SM of eq. (22) in [1]
    % Z --> F
    % U --> Y
    SM = blkvar;
    SM(1,1) = -P + Q;
    SM(1,3) = obj.A'*P - obj.C'*Y';
    SM(2,2) = -gam(k)*gam(k)*I;
    SM(2,3) = [P*F, -Y*E]';
    SM(3,3) = -P;
    SM = sdpvar(SM);

    % LMI problem to be solved
    pblmi =  [P>=0 , SM<=0]; 

    % Solve LMI conditions
    solpb = optimize(pblmi,gam(k),options_sdp);

    % Check if LMI is feasible
    if solpb.problem ~= 1
        % extract values from YALMIP symbolic decision variables
        gamma  = gam(k);
        P = value(P);
        Y = value(Y);
        break
    end
end
end


%------------- END OF CODE --------------