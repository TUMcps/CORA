function E = plusEllipsoid(E,E_cell,L,mode)
% intersectHalfspace - Computes the inner or outer approximation of the
% intersection between an ellipsoid and a halfspace
%
% Syntax:  
%    E = intersectHalfspace(E,h)
%    E = intersectHalfspace(E,h,mode)
%
% Inputs:
%    E      - ellipsoid object
%    E_cell - ellipsoid array
%    L      - directions 
%    mode   - mode ('i':inner approx; 'o': outer approx)
%
% Outputs:
%    E - ellipsoid after Minkowski sum
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      15-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
E_cell = [{E};E_cell];
N = length(E_cell);
n = dim(E);
% generally, all ellipsoids can be degenerate as long as a positively
% weighted sum of their shape matrices is invertible (otherwise all E_c are
% degenerate and and(...,'o') will not work)
% Test: If sum of Qi's is positive definite, any sum with non-zero
% weighting is (since Qi's are psd)
Q_sum = zeros(n);
TOL = min(cellfun(@(e)e.TOL,E_cell));
for i=1:N
    Q_sum = Q_sum + E_cell{i}.Q;
end

T = eye(n);
[V,S,~] = svd(Q_sum);
nt = n;
n_d = 0;
xt_rem = zeros(n_d,1);
if any(rcond(S)<=TOL)
    % all elements in E_cell are degenerate and diagonalizing with V
    % exposes common "0" dimension
    % => transform all ellipsoids so that at least one fulldimensional is
    % contained in E_cell
    T = V;
    n_d = sum(rcond(S)<=TOL);
    nt = n - n_d;
    xt_rem = zeros(n_d,1);
    for i=1:N
        Et_i = T'*E_cell{i};
        xt_rem = xt_rem + Et_i.q(nt+1:end);
        E_cell{i} = project(Et_i,1:nt);
    end
    if n_d==n
        % all Qi's zero
        E = T*ellipsoid(zeros(n),xt_rem);
        return;
    end
end    

%compute N random directions
if isempty(L)
    if nt>=2
        L = eq_point_set(nt-1,2*nt);
    else
        L = [-1,1];
    end
else
    % make sure L is normalized
    L = 1./sqrt(sum(L.^2,1)).*L;
end

E_c = lplus(E_cell,L,mode);
if strcmp(mode,'o')
    Et = and(E_c{1},E_c(2:end),mode);
else
    Et = or(E_c{1},E_c(2:end),mode);
end

% backtransform
E = T'*ellipsoid([Et.Q,zeros(nt,n_d);zeros(n_d,n)],[Et.q;xt_rem]);
%------------- END OF CODE --------------