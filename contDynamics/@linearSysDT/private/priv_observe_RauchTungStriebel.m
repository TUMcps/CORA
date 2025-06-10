function [R,tcomp] = priv_observe_RauchTungStriebel(linsysDT,params,options)
% priv_observe_RauchTungStriebel - estimates the states of linear system using
%    the Rauch-Tung-Striebel filter. This is a smoothing algorithm, which
%    cannot be used online becuase of the forward and backward propagation
%    to obtain better results than for classical state estimation. This
%    smoother is mainly intended to reconstruct states from measurements
%    as a tool for conformance checking. Unless most other implementations
%    in CORA, this observer uses concrete values and is not set-based. The 
%    implementation is close to [1]. 
%
% Syntax:
%    [R,tcomp] = priv_observe_RauchTungStriebel(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] J. M. Giron-Sierra, Digital Signal Processing with Matlab Examples,
%        Volume 3 (2017), Springer.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Stefan Liu
% Written:       15-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x0 = linsysDT.initialState; 

% Assume a process noise from set W (gaussian)
assert(all(center(params.W) == zeros(n,1)),'Center of disturbance W must be at origin')
E = generators(params.W);
Sw = 0.33^2*(E*E)'; % assume 5-sigma deviation

% Assume a process noise from set V (uniform)
assert(all(center(params.V) == zeros(q,1)),'Center of disturbance V must be at origin')
F = generators(params.V);
Sv = 0.33^2*(F*F)'; % assume 5-sigma deviation

% Kalman filter simulation preparation
%space for matrices
M = zeros(n,n); 
P = zeros(n,n,timeSteps);

%space for recording xa(n), xe(n)
rxa = zeros(n,timeSteps);
rxe = zeros(n,timeSteps);
rxa(:,1) = x0;
rxe(:,1) = x0;
xe = x0;

timerVal = tic;

%% Forward Kalman filtering
for k = 2:timeSteps
    % Prediction
    xa = linsysDT.A*xe + B*params.uTransVec(:,k-1); %a priori state
    M = linsysDT.A*P(:,:,k-1)*linsysDT.A' + Sw;
    % Update
    K = (M*linsysDT.C')*pinv(linsysDT.C*M*linsysDT.C' + Sv);
    P(:,:,k) = M - K*linsysDT.C*M;
    xe = xa + (K*(params.y(:,k)-(linsysDT.C*xa + D*params.uTransVec(:,k)))); %estimated (a posteriori) state
    % recording xa(n), xe(n)
    rxa(:,k) = xa;
    rxe(:,k) = xe;
end

%% Backward smoothing
xs = zeros(n,timeSteps);
xs(:,timeSteps) = xe; %final estimated state
M_inv = pinv(M);
for k = (timeSteps-1):-1:1
%     M = (obj.A*P(:,:,k)*obj.A')+ Sw;
    G = P(:,:,k)*linsysDT.A'*M_inv;
    xs(:,k) = rxe(:,k) + G*(xs(:,k+1)-rxa(:,k+1));
    
    % Store result
    R{k} = xs(:,k);
end

tcomp = toc(timerVal);

% ------------------------------ END OF CODE ------------------------------
