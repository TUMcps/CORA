function E = priv_andPolytope(E,P,mode)
% priv_andPolytope - computes an inner approximation or outer approximation
%    of the intersection between an ellipsoid and a polytope
%
% Syntax:
%    E = priv_andPolytope(E,P,mode)
%
% Inputs:
%    E - ellipsoid object
%    P - polytope object
%    mode - approximation of the intersection
%               'inner': inner approximation of the intersection
%               'outer': outer approximation of the intersection
%
% Outputs:
%    E - ellipsoid approximating the intersection
%
% References:
%    [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%        toolbox (ET). In Proceedings of the 45th IEEE Conference o
%        Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/and_

% Authors:       Victor Gassmann
% Written:       07-June-2022
% Last update:   05-July-2022 (VG, removed input checks; now in parent function)
% Last revision: 23-September-2024 (MW, integrate andHalfspace)

% ------------------------------ BEGIN CODE -------------------------------

% compute H-rep of P (is computed if not there; expensive!)
[A,b] = constraints(P);

% loop over each inequality constraint
for i=1:size(A,1)
    E = aux_andHalfspace(E,polytope(A(i,:),b(i)),mode);
end

end


% Auxiliary functions -----------------------------------------------------

function E = aux_andHalfspace(E,P,mode)

% compute distance to corresponding hyperplane
dist = distance(E,polytope([],[],P.A,P.b));
n = dim(E);

% touching, completely inside or outside
if dist >= -E.TOL
    if P.A*E.q > P.b
        % completely outside or touching
        if withinTol(dist,0,E.TOL)
            % touching
            % point on hyperplane
            xh = P.A'*P.b;
            % direction resulting in touching point has positive 
            % inner product with xh
            v = sign((xh-E.q)'*P.A')*P.A';
            % compute touching point
            [~,x] = supportFunc_(E,v,'upper');
            E = ellipsoid(zeros(n),x);
        else 
            E = ellipsoid.empty(n);
        end
    else
        % E completely inside (or touching)
        % E = E;
    end
    return;
end
% ...now established that they are intersecting

T = eye(n);
x_rem = [];

if ~isFullDim(E)
    nt = rank(E);
    % check if E.Q all zero
    if nt==0
        % check if E is contained in the polytope
        if contains_(P,E.q,'exact',0)
            E = ellipsoid(zeros(n),E.q);
        else
            E = ellipsoid.empty(n);
        end
        return;
    end
    [T,~,~] = svd(E.Q);
    E = T'*E;
    % transform inequality constraint (possible since T invertible)
    P = polytope(P.A*T,P.b);
    % project
    x_rem = E.q(nt+1:end);
    E = project(E,1:nt);
    P = polytope(P.A(1:nt),P.b-P.A(nt+1:end)*x_rem);
end

n_nd = dim(E);
% normalize inequality constraint
% shift E and P such that P.b = 0 and transform such that c=e_1
A = P.A'/norm(P.A); b = P.b/norm(P.A);
% compute transformation matrix so that e_1 = S*c;
unit_vector_1 = unitvector(1,n_nd);
S = vecalign(unit_vector_1,A);
P = polytope(unit_vector_1',0);
E = -b*unit_vector_1 + S*E;

n_rem = n-n_nd;
% now non-degenerate
W1 = inv(E.Q);
q1 = E.q;


if strcmp(mode,'outer')
    [r_s,~] = supportFunc_(E,unit_vector_1,'lower');
    % makes more sense than ET original: define degenerate ellipsoid that
    % covers the transformed ellipsoid "exactly"
    q2 = [1/2*r_s;zeros(n_nd-1,1)];
    W2 = diag([4/r_s^2;zeros(n_nd-1,1)]);
    % also, ET original does not work?
    
    p = priv_compIntersectionParam(W1,q1,W2,q2);
    [~,Q_nd,q_nd] = priv_rootfnc(p,W1,q1,W2,q2);
else
    % that is "ellipsoidal toolbox original" (not sure why this works)
    E_hyp = and_(E,P,'outer');
    q2 = E_hyp.q-2*sqrt(max(eig(E.Q)))*unit_vector_1;
    W2 = (unit_vector_1*unit_vector_1')*1/(4*max(eig(E.Q)));
    b1 = (E.q-E_hyp.q)'*W1*(E.q-E_hyp.q);
    [~,xb] = supportFunc_(E,-unit_vector_1,'upper');
    b2 = (q2-xb)'*W2*(q2-xb);
    %
    b1 = min(1,b1);
    b2 = min(1,b2);
    assert(b2<1,'Cannot be =1, since then W is not invertible');
    t1 = (1-b2)/(1-b1*b2);
    t2 = (1-b1)/(1-b1*b2);
    q_nd = (t1*W1+t2*W2)\(t1*W1*q1+t2*W2*q2);
    W = t1*W1+t2*W2;
    Q_nd = (1-t1*q1'*W1*q1-t2*q2'*W2*q2+q_nd'*W*q_nd)*inv(W);
end
E_nd = ellipsoid(Q_nd,q_nd);
% revert S transform + shift
Et = S'*E_nd + dist*A;

% restore original dimensions and backtransform
E_t = ellipsoid([Et.Q,zeros(n_nd,n_rem);zeros(n_rem,n)],[Et.q;x_rem]);
E = T*E_t;

end

% ------------------------------ END OF CODE ------------------------------
