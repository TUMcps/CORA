function [R,tcomp] = observe_RauchTungStriebel(obj,options)
% observe_RauchTungStriebel - estimates the states of linear
% system using the Rauch-Tung-Striebel filter. This is a smoothing 
% algorithm, which cannot be used online becuase of the forward and backward
% propagation to obtain better results than for classical state estimation.
% This smoother is mainly intended to reconstruct states from measurements
% as a tool for conformance checking. Unless most other implementations in 
% CORA, this observer uses concrete values and is not set-based. The 
% implementation is close to [1]. 
%
% Syntax:
%    [R,tcomp] = observe_RauchTungStriebel(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
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
% Example:
%    -
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

x0 = obj.initialState; 

% Assume a process noise from set W (gaussian)
assert(all(W.center == zeros(n,1)),'Center of disturbance W must be at origin')
E = W.generators;
Sw = 0.33^2*(E*E)'; % assume 5-sigma deviation

% Assume a process noise from set V (uniform)
assert(all(V.center == zeros(q,1)),'Center of disturbance V must be at origin')
F = V.generators;
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

tic

%% Forward Kalman filtering
for k = 2:timeSteps
    % Prediction
    xa = obj.A*xe + B*options.uTransVec(:,k-1); %a priori state
    M = obj.A*P(:,:,k-1)*obj.A' + Sw;
    % Update
    K = (M*obj.C')*pinv(obj.C*M*obj.C' + Sv);
    P(:,:,k) = M - K*obj.C*M;
    xe = xa + (K*(options.y(:,k)-(obj.C*xa + D*options.uTransVec(:,k)))); %estimated (a posteriori) state
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
    G = P(:,:,k)*obj.A'*M_inv;
    xs(:,k) = rxe(:,k) + G*(xs(:,k+1)-rxa(:,k+1));
    
    % Store result
    R{k} = xs(:,k);
end

tcomp = toc;

end


% ------------------------------ END OF CODE ------------------------------
