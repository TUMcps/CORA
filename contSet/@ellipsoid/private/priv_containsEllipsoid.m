function [res,cert,scaling] = priv_containsEllipsoid(E1,E2,tol)
% priv_containsEllipsoid - checks whether an ellipsoid contains another
%    ellipsoid
%
% Syntax:
%    E = priv_containsEllipsoid(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object (circumbody)
%    E2 - ellipsoid object (inbody)
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, E2 is
%           guaranteed to not be contained in E1, whereas if res=false and
%           cert=false, nothing can be deduced (E2 could still be
%           contained in E1).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(E1 - E1.center) + E1.center contains E2.
%
% References:
%    [1] Boyd et al. Convex Optimization (B.2, ex. B.1)
%    [2] Kulmburg A., Schäfer L., Althoff A., Approximability of the
%          Containment Problem for Zonotopes and Ellipsotopes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/contains_

% Authors:       Victor Gassmann, Adrian Kulmburg
% Written:       16-March-2021
% Last update:   27-July-2021 (allowed both E1,E2 to be degenerate)
%                23-May-2022 (solver-specific implementation of SDP to avoid yalmip)
%                06-July-2022 (VG, rescaling of M/Ms to avoid numerical issues)
%                20-January-2025 (AK, overhauled the algorithm, aligning it with [2]) 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% We follow the naming conventions from [2]
n = dim(E1);
G = generators(E2);
H = generators(E1);

% Need to check if E1 even has a chance to contain E2, in case E1 is
% degenerate
r = rank(H, tol);
if r < n
    r_comb = rank([H G], tol);
    if r_comb > r
        res = false;
        cert = true;
        scaling = inf;
        return
    end
end
% In any other instance, we should be good

m = size(G, 2);
ell = size(H, 2);

Theta = pinv(H) * G;
theta = pinv(H) * (E2.q - E1.q);

rho = sdpvar(1,1,'full');
delta = sdpvar(1,1,'full');

cost = rho;
A = [Theta'*Theta Theta'*theta; theta'*Theta theta'*theta];
B = zeros([m+1 m+1]);
B(m+1,m+1) = 1;
C = eye(m+1);
C(m+1,m+1) = -1;
constraint = [A-rho*B<=delta*C  delta>=0];

% solve with solver
if isSolverInstalled('mosek')
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
else
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
end
diagnostics = optimize(constraint, cost, options);

rho = value(rho);
if rho <= 1+tol
    res = true;
else
    res = false;
end
cert = true;
scaling = sqrt(rho);

% ------------------------------ END OF CODE ------------------------------
