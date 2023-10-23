function val = distanceEllipsoid(E1,E2)
% distanceEllipsoid - computes the distance between two ellipsoids
%
% Syntax:
%    val = distanceEllipsoid(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object
%
% Outputs:
%    val - distance between ellipsoids
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   19-May-2022
%                02-June-2022 (VG, complete rework)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% not supported for degenerate ellipsoids or of different size
if dim(E1) ~= dim(E2)
    throw(CORAerror('CORA:dimensionMismatch',E1,E2));
end

if ~isFullDim(E1) && ~isFullDim(E2)
    throw(CORAerror('CORA:degenerateSet',...
        'At least one ellipsoid has to be non-degenerate!'));
end

% shift E1 to zero (not strictly necessary)
q1 = E1.q;
E1 = E1 + (-q1);
E2 = E2 + (-q1);

x2_rem = [];
% make E1 always non-degenerate
if ~isFullDim(E1)
    tmp = E1;
    E1 = E2;
    E2 = tmp;
end
n = dim(E1);
nt = n;

if ~isFullDim(E2)
    nt = rank(E2);
    % Q all zero
    if nt==0
        % distance from E1 to center q
        val = distancePoint(E1,E2.q);
        return;
    end
    [T,~,~] = svd(E2.Q);
    % transform E1 and E2 such that degeneracy of E2 is axis-aligned
    E2 = T'*E2;
    E1 = T'*E1;
    x2_rem = E2.q(nt+1:end);
    % project
    E2 = project(E2,1:nt);
end

% solve 
    %%% min_{x}   (x(1:nt)-E2.q)'*inv(E2.Q)*(x(1:nt)-E2.q) + 
    %%%           (x(nt+1:n)-x2_rem)'*1/TOL*(x(nt+1:n)-x2_rem),
    %%%     s.t.  (x-E1.q)'*inv(E1.Q)*(x-E1.q) <= 1.

if isSolverInstalled('coneprog')
    % reformulate to (equivalently, i.e. objval is sqrt of wanted value)
    %%% min_{z,t1,t2,x} z,
    %%%     s.t.        ||[t1;t2]||_2 <= z,
    %%%                 ||sqrtm(E2.Q)\(x(1:nt)-E2.q)||_2 <= t1,
    %%%                 ||1/sqrt(E1.TOL)*(x(nt+1:n)-x2_rem)||_2 <= t2,
    %%%                 ||sqrtm(E1.Q)\(x-E1.q)||_2 <= 1.

    % opt var = [z;t1;t2;x]
    % cone for z
    A_z = [zeros(2,1),eye(2),zeros(2,n)];
    b_z = zeros(size(A_z,1),1);
    d_z = [1;zeros(1+1+n,1)];
    gamma_z = 0;
    socs(4) = secondordercone(A_z,b_z,d_z,gamma_z);

    % cone for t1
    B2_inv = inv(sqrtm(E2.Q));
    A_t1 = [zeros(nt,1+1+1),B2_inv,zeros(nt,n-nt)];
    b_t1 = B2_inv*E2.q;
    d_t1 = [0;1;0;zeros(n,1)];
    gamma_t1 = 0;
    socs(3) = secondordercone(A_t1,b_t1,d_t1,gamma_t1);

    % cone for t2
    A_t2 = [zeros(n-nt,1+1+1),zeros(n-nt,nt),1/sqrt(E1.TOL)*eye(n-nt)];
    b_t2 = 1/sqrt(E1.TOL)*x2_rem;
    d_t2 = [0;0;1;zeros(n,1)];
    gamma_t2 = 0;
    socs(2) = secondordercone(A_t2,b_t2,d_t2,gamma_t2);

    % cone for E1
    B1_inv = inv(sqrtm(E1.Q));
    A_E1 = [zeros(n,1+1+1),B1_inv];
    b_E1 = B1_inv*E1.q;
    d_E1 = zeros(1+1+1+n,1);
    gamma_E1 = -1;

    socs(1) = secondordercone(A_E1,b_E1,d_E1,gamma_E1);
    
     % supress output
    options = optimoptions(@coneprog,'display','none');

    % solve problem
    [~,objval_] = coneprog([1;zeros(1+1+n,1)],socs,[],[],[],[],[],[],options);
    objval = objval_^2;
    val = objval - 1;

else
    % solve 
    %%% min_{x}   (x(1:nt)-E2.q)'*inv(E2.Q)*(x(1:nt)-E2.q) + 
    %%%           (x(nt+1:n)-x2_rem)'*1/TOL*(x(nt+1:n)-x2_rem),
    %%%     s.t.  (x-E1.q)'*inv(E1.Q)*(x-E1.q) <= 1.
    
    fobj = @(x) norm(sqrtm(E2.Q)\(x(1:nt)-E2.q))^2 + ...
                norm(1/sqrt(E1.TOL)*(x(nt+1:n)-x2_rem))^2;
    gobj = @(x) 2*[E2.Q\(x(1:nt)-E2.q);1/E1.TOL*(x(nt+1:n)-x2_rem)];
    f_obj = @(x) deal(fobj(x),gobj(x));
    fcstr = @(x) norm(sqrtm(E1.Q)\(x-E1.q))^2 - 1;
    gcstr = @(x) 2*(E1.Q\(x-E1.q));
    f_cstr = @(x) deal(fcstr(x),[],gcstr(x),[]);
    H = @(x,lambda) 2*blkdiag(inv(E2.Q),1/E1.TOL*eye(n-nt)) + ...
                        lambda.ineqnonlin(1)*2*inv(E1.Q);
    
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
                            'display','none','SpecifyObjectiveGradient',true,...
                            'SpecifyConstraintGradient',true,...
                            'HessianFcn',H);
    
    [~,objval] = fmincon(f_obj,E1.q,[],[],[],[],[],[],f_cstr,options);
    val = objval - 1;
end

% ------------------------------ END OF CODE ------------------------------
