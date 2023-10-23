function [OGain,P,gamma,lambda,tComp] = observe_gain_ESO_C(obj,options)
% observe_gain_ESO_C - computes the gain for the guaranteed state estimation
% approach from [1]. 
%
%
% Syntax:
%    [OGain,tComp]= observe_gain_ESO_C(obj,options)
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
%    [1] Nassim Loukkas, John J. Martinez, and Nacim Meslem. Set-
%        membership observer design based on ellipsoidal invariant
%        sets. IFAC-Papers On Line, 50(1):6471-6476, 2017.
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
% Last update:   16-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic

% obtain system dimension and nr of outputs
n = obj.dim; 

% set options of solver
options_sdp = sdpsettings;
options_sdp.solver = options.solver;
options_sdp.verbose = 1;  

%% Alg. 1 of [1] 
% initialization of Alg. 1
Q = eye(n); 

% line 2-6 of Alg. 1 in [1]
[~, ~, ~, ~, Ao] = aux_compute_Q(obj,Q,options,options_sdp);

% Compute the disturbance variance W; we use a more precise method compared
% to using (27) in [1]. Quantile function for probability p of the 
% chi-squared distribution so that the probability of being in the
% confidence set is p
p = 0.999;
quantileFctValue = chi2inv(p,dim);

% obtain covariance matrix
W = options.W.Q/quantileFctValue;

% Obtain the covariance matrix V using (26) in [1]
V = sdpvar(n, n, 'symmetric');
Fa = [V>=0, Ao*V*Ao'-V<=-W];    
optimize(Fa);    
V = value(V);

% obtain Q = V^-1
Q = inv(V);

% go back to line 2-6 of Alg. 1 in [1]
[Q, P, Y, gamma] = aux_compute_Q(obj,Q,options,options_sdp);

% Compute λ as the minimum generalized eigenvalue of the pair (Q,P), line 7
lambda = min(eigs(Q,P));

% store gain
OGain = P\Y;

% computation time
tComp = toc;

end


% Auxiliary functions -----------------------------------------------------

% This function realizes line 2-6 of Alg. 1 in [1]
function [Q, P, Y, gamma, Ao] = aux_compute_Q(obj,Q,options,options_sdp)
    % solve LMI, line 2 of Alg. 1
    [P,Y,gamma] = aux_solveLMI(obj,Q,options,options_sdp);

    % compute L, line 3 of Alg. 1
    L = P\Y;

    % Compute observer matrix, line 4 of Alg. 1
    Ao = obj.A-L*obj.C;

    % Compute E, line 5 of Alg. 1
    E = [options.W.Q, -L*options.V.Q];

    %% Using P and γ, find a new matrix Q which satisfies the LMI in (20) 
    % of [1] and minimizes (− log(det(Q)))
    
    % obtain system dimension and nr of outputs
    n = obj.dim; 
    nrOfOutputs = size(obj.C,1);
    nrOfDistGens = size(options.W.Q,2);
    
    % identity matrix
    I = eye(nrOfDistGens + nrOfOutputs);
    
    % create symmetric matrix SM of LMI problem
    Q = sdpvar(n,n,'symmetric');
    SM = blkvar;
    SM(1,1) = Ao'*P*Ao-P+Q;
    SM(1,2) = Ao'*P*E;
    SM(2,1) = E'*P*obj.A;
    SM(2,2) = E'*P*E- gamma^2*I;
    SM = sdpvar(SM);
    constraint = [Q>=0, SM<=0]; % constraint   
    objective = -log(det(Q));   % objective function
    ops_lmi = sdpsettings('solver','penlab','verbose',0,'warning',0);    
    solpb = optimize(constraint, objective, ops_lmi); % optimze the LMIs
    % return Q
    Q = double(Q);
end


% Auxiliary functions -----------------------------------------------------

function [P,Y,gamma] = aux_solveLMI(obj,Q,options,options_sdp)
    
    % shape matrices of the disturbance and porcess noise sets
    F = options.W.Q;
    E = options.V.Q;
    
    % obtain system dimension and nr of outputs
    n = obj.dim;
    nrOfOutputs = size(obj.C,1);
    nrOfDistGens = size(F,2);

    %% define YALMIPs symbolic decision variables
    % state
    P = sdpvar(n,n,'symmetric'); 
    % gain matrix
    Y = sdpvar(n,nrOfOutputs,'full'); 
    % identity matrix
    I = eye(nrOfDistGens + nrOfOutputs);

    % possible values of gam
    gam = linspace(1,50,100);

    % Optimization loop
    for k = 1:length(gam)
        % implementation of symmetric matrix SM of eq. (19) in [1]
        % Z --> E
        % U --> Y
        SM = blkvar;
        SM(1,1) = -P + Q;
        SM(1,3) = obj.A'*P - obj.C'*Y';
        SM(2,2) = -gam(k)*gam(k)*I;
        SM(2,3) = ([P*F, - Y*E])';
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


% ------------------------------ END OF CODE ------------------------------
