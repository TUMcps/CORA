function val = distanceMptPolytope(E,P)
% distanceMptPolytope - computes the distance from an ellipsoid to a
% mptPolytope object
%
% Syntax:  
%    val = distanceMptPolytope(E,H)
%
% Inputs:
%    E - ellipsoid object
%    P - mptPolytope object
%
% Outputs:
%    val - distance between E and P
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  18-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
n = dim(E);
x_rem = zeros(0,1);
% check if ellipsoid is degenerate
if E.isdegenerate
    nt = rank(E);
    % check if E.Q is all zero
    if nt==0
        val = distance(P,E.q);
        return;
    end
    [T,~,~] = svd(E.Q);
    E = T'*E;
    % transform mptPolytope
    P = T'*P;
    x_rem = E.q(nt+1:end);
    % project
    E = project(E,1:nt);
end
n_nd = dim(E);

% normalize to prevent issues
A = P.P.A;
b = P.P.b;
% normalize halfspace rep
fac = 1./sqrt(sum(A.^2,2));
A = fac.*A;
b = fac.*b;

% solve optimization problem (solution is >=0)
x_nd = sdpvar(n_nd,1);
x = [x_nd;x_rem];
y = sdpvar(n,1);
C_e = (x_nd-E.q)'*inv(E.Q)*(x_nd-E.q)<=1;
C_p = A*y<=b;
f_obj = norm(x-y);
optimize([C_e,C_p],f_obj,sdpsettings('verbose',0));
% extract distance (is the same as in untransformed space since T is
% unitary)
val = value(f_obj);
%------------- END OF CODE --------------