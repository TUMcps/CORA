function val = distancePolytope(E,P)
% distancePolytope - computes the distance from an ellipsoid to a
%    polytope
%
% Syntax:
%    val = distancePolytope(E,H)
%
% Inputs:
%    E - ellipsoid object
%    P - polytope object
%
% Outputs:
%    val - distance between ellipsoid and polytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   18-March-2021
%                03-June-2022 (VG, complete rework)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(E);
x_rem = zeros(0,1);
nt = n;
% check if ellipsoid is degenerate
if ~isFullDim(E)
    nt = rank(E);
    % check if E.Q is all zero
    if nt==0
        val = distance(P,E.q);
        return;
    end
    [T,~,~] = svd(E.Q);
    E = T'*E;
    % transform polytope
    P = T'*P;
    x_rem = E.q(nt+1:end);
    % project
    E = project(E,1:nt);
end

% normalize to prevent issues
A = P.A;
b = P.b;
% normalize halfspace rep
fac = 1./sqrt(sum(A.^2,2));
A = fac.*A;
b = fac.*b;

% solve 
%%% min_{x}   (x(1:nt)-E.q)'*inv(E.Q)*(x(1:nt)-E.q) + 
%%%           (x(nt+1:n)-x_rem)'*1/TOL*(x(nt+1:n)-x_rem),
%%%     s.t.  A*x <= b.

% solve using quadprog
H = 2*blkdiag(inv(E.Q),1/E.TOL*eye(n-nt));
f = -2*[E.Q\E.q;1/E.TOL*x_rem];

% supress output
options = optimoptions(@quadprog,'display','none');

[~,objval_] = quadprog(H,f,A,b,[],[],[],[],[],options);

objval = objval_ + norm(sqrtm(E.Q)\E.q)^2 + norm(1/sqrt(E.TOL)*x_rem)^2;
val = objval - 1;

% if isSolverInstalled('coneprog')
%     % reformulate to (equivalently, i.e. objval is sqrt of wanted value)
%     %%% min_{z,t1,t2,x} z,
%     %%%     s.t.        ||[t1;t2]||_2 <= z,
%     %%%                 ||sqrtm(E.Q)\(x(1:nt)-E.q)||_2 <= t1,
%     %%%                 ||1/sqrt(E1.TOL)*(x(nt+1:n)-x_rem)||_2 <= t2,
%     %%%                 A*x <= b.
% 
%     % opt var = [z;t1;t2;x]
%     % cone for z
%     A_z = [zeros(2,1),eye(2),zeros(2,n)];
%     b_z = zeros(size(A_z,1),1);
%     d_z = [1;zeros(1+1+n,1)];
%     gamma_z = 0;
%     socs(3) = secondordercone(A_z,b_z,d_z,gamma_z);
% 
%     % cone for t1
%     B_inv = inv(sqrtm(E.Q));
%     A_t1 = [zeros(nt,1+1+1),B_inv,zeros(nt,n-nt)];
%     b_t1 = B_inv*E.q;
%     d_t1 = [0;1;0;zeros(n,1)];
%     gamma_t1 = 0;
%     socs(2) = secondordercone(A_t1,b_t1,d_t1,gamma_t1);
% 
%     % cone for t2
%     A_t2 = [zeros(n-nt,1+1+1),zeros(n-nt,nt),1/sqrt(E.TOL)*eye(n-nt)];
%     b_t2 = 1/sqrt(E.TOL)*x_rem;
%     d_t2 = [0;0;1;zeros(n,1)];
%     gamma_t2 = 0;
%     socs(1) = secondordercone(A_t2,b_t2,d_t2,gamma_t2);
%     
%     % extend linear inequalities to new opt vars
%     A = [zeros(size(A,1),1+1+1),A];
% 
%      % supress output
%     options = optimoptions(@coneprog,'display','none');
% 
%     % solve problem
%     [~,objval_] = coneprog([1;zeros(1+1+n,1)],socs,A,b,[],[],[],[],options);
%     objval = objval_^2;
%     val = objval - 1;
% end

% ------------------------------ END OF CODE ------------------------------
