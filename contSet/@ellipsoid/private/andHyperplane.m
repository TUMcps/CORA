function E = andHyperplane(E,H)
% andHyperplane - computes the exact intersection of an ellipsoid and a
%    hyperplane
%
% Syntax:
%    E = andHyperplane(E,H)
%
% Inputs:
%    E - Ellipsoid object
%    H - conHyperplane object (but {x|C*x<=d}= R^n)
%
% Outputs:
%    E - ellipsoid object
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if distance(E,H)>E.TOL
    E = ellipsoid;
    return;
end

% they are intersecting
n = length(E.q);
x_rem = zeros(0,1);
T = eye(n);

if ~isFullDim(E)
    nt = rank(E);
    % check if E.Q is all zero
    if nt==0
        % if E is 0-d, the result is either E.q if E.q\in H, or empty set
        if contains_(H,E.q,'exact',0)
            E = ellipsoid(zeros(n),E.q);
        else
            E = ellipsoid;
        end
        return;
    end
    [T,~,~] = svd(E.Q);
    % transform E such that degenerate dimensions are axis-aligned
    E = T'*E;
    % transform hyperplane (possible since T unitary)
    H = conHyperplane((T'*H.a')',H.b);
    % project
    x_rem = E.q(nt+1:end);
    E = project(E,1:nt);
    H = conHyperplane(H.a(1:nt),H.b-H.a(nt+1:end)*x_rem);
end

n_nd = length(E.q);
n_rem = n-n_nd;
% normalize hyperplane
c = H.a'/norm(H.a');
d = H.b/norm(H.a');
H = conHyperplane(c',d);
I = eye(n_nd);

% check if only 1 non-degenerate dimension remains in E_nd
if n_nd==1
    % E, H are 1d: check if intervals intersect
    % compute enclosing interval
    IntE = E.q + interval(-sqrt(E.Q),sqrt(E.Q));
    xH = H.b/H.a';
    x_max = max(abs(xH));
    r_xH = x_max*E.TOL;
    IntE_TOL = IntE + interval(-r_xH,r_xH);
    if ~contains_(IntE_TOL,xH,'exact',0)
        E = ellipsoid;
        return;
    end
    E_t = ellipsoid(0,xH);
else
    % compute transformation matrix so that I(:,1) = S*c;
    
    % for more detail on the following, see [1]
    S = vecalign(I(:,1),c);
    E = -I(:,1)*d+S*E;
    % transformed hyperplane
    H = conHyperplane(I(:,1)',0);
    M = inv(E.Q);
    Mb = M(2:end,2:end);
    mb = M(2:end,1);
    m11 = M(1,1);
    Mbinv = inv(Mb);
    w_s = E.q+E.q(1)*[-1;Mbinv*mb];
    a = 1-E.q(1)^2*(m11-mb'*Mbinv*mb);
    if a<0 && a>-E.TOL
        a = 0;
    elseif a<-E.TOL
        throw(CORAerror('CORA:specialError',...
            'Error computing intersection of ellipsoid and hyperplane!'));
    end
    W_s = a*[zeros(1,n_nd);[zeros(n_nd-1,1),Mbinv]];
    Ew = ellipsoid(W_s,w_s);
    E_t = S'*Ew+d*c;
end

% reintroduce x_rem and backtransform
E = ellipsoid([E_t.Q,zeros(n_nd,n_rem);zeros(n_rem,n)],[E_t.q;x_rem]);
E = T*E;

% ------------------------------ END OF CODE ------------------------------
