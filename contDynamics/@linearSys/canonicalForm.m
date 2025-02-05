function [linsys_,U_,u_,V_,v_] = canonicalForm(linsys,U,uVec,W,V,vVec)
% canonicalForm - rewrite inhomogeneity to canonical forms:
%    Ax + Bu + c + w  ->  Ax + u', where u' \in U_ + u_
%    Cx + Du + k + v  ->  Cx + v', where v' \in V_ + v_
%
% Syntax:
%    [linsys_,U_,u_,V_,v_] = canonicalForm(linsys,U,uVec,W,V,vVec)
%
% Inputs:
%    linsys - linearSys object
%    U - input set (time-varying)
%    uVec - input vector (piecewise constant)
%    W - disturbance set (time-varying)
%    V - noise set (time-varying)
%    vVec - noise vector (piecewise constant)
%
% Outputs:
%    sys_ - linearSys object in canonical form
%    U_ - input set, interpreted as time-varying over t
%    u_ - piecewise constant input trajectory
%    V_ - output uncertainty
%    v_ - piecewise constant offset on output
%
% Example:
%    linsys = linearSys([-1 -4; 4 -1],[1; -1]);
%    U = zonotope(1,0.02);
%    uVec = [2 4 -1 0 3];
%    W = zonotope(0.5,0.05);
%    V = zonotope(-0.2,0.1);
%    vVec = [0 -1 0 1 2 3];
%
%    [linsys_,U_,u_,V_,v_] = canonicalForm(linsys,U,uVec,W,V,vVec);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/linearSys

% Authors:       Mark Wetzlinger
% Written:       05-April-2024
% Last update:   17-October-2024 (MW, integrate disturbance/noise matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% given sets:
% - U: set, not necessarily centered at zero
% - uVec: vector or list of vectors
% - W: set, not necessarily centered at zero
% - V: set, not necessarily centered at zero
% - vVec: vector or list of vectors

% read out all centers
centerU = center(U);
centerW = center(W);
centerV = center(V);

% shift all sets so that they are centered at the origin
U = U - centerU;
W = W - centerW;
V = V - centerV;

% offset vector the output: if it is not constant, we require an additional
% column in U because there are steps+1 time points (on which we evaluate
% the output set) but only steps time intervals (for which we need uVec)
if size(uVec,2) == 1
    v_ = linsys.D * uVec + linsys.k + linsys.F * (centerV + vVec);
elseif ~any(linsys.D,'all') && ~any(linsys.k) && ~any(linsys.F,'all')
    v_ = zeros(linsys.nrOfOutputs,1);
else
    % only compute if a non-zero result is to be expected
    v_ = linsys.D * [uVec, zeros(linsys.nrOfInputs,1)] + linsys.k + linsys.F * (centerV + vVec);
end
% simplify representation if result is all-zero
if ~any(v_,'all')
    v_ = zeros(linsys.nrOfOutputs,1);
end
% time-varying uncertainty on the output
V_ = linsys.D * (U + centerU) + linsys.F * V;

% offset vector for state
u_ = linsys.B * uVec + linsys.B * centerU + linsys.c + linsys.E * centerW;
% time-varying uncertainty for state
U_ = linsys.B * U + linsys.E * W;


% update system dynamics
n = linsys.nrOfDims; r = linsys.nrOfOutputs;
linsys_ = linearSys(linsys.A,eye(n),zeros(n,1),...       % A, B, c
                linsys.C,zeros(r,n),zeros(r,1),...   % C, D, k
                zeros(n),eye(r));                 % E, F
% copy helper property
linsys_.taylor = linsys.taylor;

% ------------------------------ END OF CODE ------------------------------
