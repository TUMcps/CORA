function [linsys,R_uTrans_apx] = priv_krylov_R_uTrans(linsys,t,uTrans,options)
% priv_krylov_R_uTrans - computes the reachable continuous set for time t 
%                   due to constant input vector uTrans
%
% Syntax:
%    [linsys,R_uTrans_apx] = priv_krylov_R_uTrans(linsys,t,uTrans,options)
%
% Inputs:
%    linsys - linearSys object modeling the integrator system for inputs
%             (see [1] Lemma 4.1) 
%    t      - time point (scalar double)
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    linsys   - original linearSys object (with results of this function saved)
%    R_uTrans_apx - reachable set due to uTrans at the time t
%
% References:
%    [1] Matthias Althoff, Reachability Analysis of Large Linear Systems
%        With Uncertain Inputs in the Krylov Subspace
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: - 

% Authors:       Maximilian Perschl
% Written:       25-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

state_dim = size(linsys.A,1);

% create integrator system
A_int = [linsys.A, uTrans; zeros(1,state_dim), 0];

% equivalent initial state
eqivState = [zeros(state_dim,1); 1];

% Arnoldi
[V_uT,H_uT,~] = ...
    priv_subspace_Krylov_individual_Jawecki(A_int,eqivState,1,options);

V_uT_proj = linsys.C*V_uT(1:state_dim,:);

%initialize linear reduced dynamics
linRedSys = linearSys('reduced_sys_uTrans',H_uT,V_uT_proj');
% need taylor object for correction matrix
linRedSys.taylor = taylorLinSys(H_uT);

% compute P*(...)*e_1 via indexing for fast computation
% norm(equiv_state) = 1, so we leave it out;
expMat = expm(H_uT*t);
R_uTrans_apx = V_uT_proj * expMat(:,1);

% add krylov error
if options.krylovError > 2*eps
    R_uTrans_proj = zonotope(R_uTrans_apx,eye(size(R_uTrans_apx,1))*options.krylovError);
else
    R_uTrans_proj = zonotope(R_uTrans_apx);
end

G = priv_correctionMatrixInput(linRedSys,t,options.taylorTerms);
% compute P*(...)*e_1 via indexing for fast computation
% norm(equiv_state) = 1, so we leave it out;
inputCorr_unprojected = V_uT(1:state_dim,:)*G(:,1);
% inputCorrection = zonotope(inputCorr_unprojected(1:state_dim)); 
inputCorr_unprojected = zonotope(inputCorr_unprojected); 


% save results for future re-computations of correction matrix
linsys.krylov.uTrans_sys = linRedSys;
linsys.krylov.R_uTrans = V_uT(1:state_dim,:)*expMat(:,1);
linsys.krylov.R_uTrans_proj = R_uTrans_proj;
linsys.krylov.inputCorr = linsys.C*inputCorr_unprojected;
linsys.krylov.inputCorr_radius = radius(inputCorr_unprojected);

end

% ------------------------------ END OF CODE ------------------------------
