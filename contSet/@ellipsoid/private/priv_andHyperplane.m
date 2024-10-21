function E = priv_andHyperplane(E,P)
% priv_andHyperplane - computes the exact intersection of an ellipsoid and
%    a hyperplane
%
% Syntax:
%    E = priv_andHyperplane(E,P)
%
% Inputs:
%    E - ellipsoid object
%    P - polytope object representing a hyperplane
%
% Outputs:
%    E - ellipsoid representing the intersection
%
% References: 
%    [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%        https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/and_

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension
n = dim(E);

% check if the intersection between the ellipsoid and the hyperplane is
% empty using distance
if distance(E,P) > E.TOL
    E = ellipsoid.empty(n);
    return;
end

% check for degeneracy
isDeg = ~isFullDim(E);
if isDeg
    n_subspace = rank(E);
    % check if E.Q is all zero
    if n_subspace==0
        % if E is 0-d, the result is either E.q if E.q\in H, or empty set
        if contains_(P,E.q,'exact',0)
            E = ellipsoid(zeros(n),E.q);
        else
            E = ellipsoid.empty(n);
        end
        return;
    end
    [T,~,~] = svd(E.Q);
    % transform E such that degenerate dimensions are axis-aligned
    E = T'*E;
    % transform hyperplane (possible since T unitary)
    P = polytope([],[],(T'*P.Ae')',P.be);
    % project
    x_rem = E.q(n_subspace+1:end);
    E = project(E,1:n_subspace);
    P = polytope([],[],P.Ae(1:n_subspace),...
                       P.be-P.Ae(n_subspace+1:end)*x_rem);

    n_rem = n-n_subspace;
end

n_nd = dim(E);
% normalize hyperplane
P_ = normalizeConstraints(P,'A');

% check if only 1 non-degenerate dimension remains in E_nd
if n_nd==1
    % ellipsoid and hyperplane are 1D: check if intervals intersect

    % compute enclosing interval
    IntE = E.q + interval(-sqrt(E.Q),sqrt(E.Q));
    xH = P_.be / P_.Ae';

    r_xH = max(abs(xH)) * E.TOL;
    IntE_TOL = IntE + interval(-r_xH,r_xH);
    if ~contains_(IntE_TOL,xH,'exact',0)
        E = ellipsoid.empty(n);
        return;
    end
    E_t = ellipsoid(0,xH);

else
    % compute transformation matrix so that e_1 = S*P.be;
    
    % for more detail on the following, see [1]
    unitvector_1 = unitvector(1,n_nd);
    S = vecalign(unitvector_1,P_.Ae');
    E = -unitvector_1*P_.be + S*E;

    % transformed hyperplane (not needed?)
    % P_ = polytope([],[],unitvector_1',0);
    M = inv(E.Q);
    Mb = M(2:end,2:end);
    mb = M(2:end,1);
    m11 = M(1,1);
    Mbinv = inv(Mb);
    w_s = E.q+E.q(1)*[-1;Mbinv*mb];
    a = 1-E.q(1)^2*(m11-mb'*Mbinv*mb);
    if a < 0 && a > -E.TOL
        a = 0;
    elseif a < -E.TOL
        throw(CORAerror('CORA:specialError',...
            'Error computing intersection of ellipsoid and hyperplane!'));
    end

    W_s = a*[zeros(1,n_nd);[zeros(n_nd-1,1),Mbinv]];
    Ew = ellipsoid(W_s,w_s);
    E_t = S'*Ew + P_.be*P_.Ae';
end

% degenerate case: reintroduce x_rem and backtransform
if isDeg
    E = ellipsoid([E_t.Q,zeros(n_nd,n_rem);zeros(n_rem,n)],[E_t.q;x_rem]);
    E = T*E;
else
    E = E_t;
end

% ------------------------------ END OF CODE ------------------------------
