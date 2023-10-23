function E = plusEllipsoid(E,L,mode)
% plusEllipsoid - Computes an inner-approximation or outer-approximation of
%    the Minkowski sum of a list of ellipsoids
%
% Syntax:
%    E = plusEllipsoid(E,L,mode)
%
% Inputs:
%    E      - ellipsoid object (or array)
%    L      - directions 
%    mode   - type of approximation: 'inner', 'outer'
%
% Outputs:
%    E - ellipsoid object after Minkowski sum
%
% References:
%   [1]: A. Halder. "On the parameterized computation of minimum volume
%        outer ellipsoid of Minkowski sum of ellipsoids.", CDC 2018. 
%   [2]: S. Boyd et al. "Linear matrix inequalities in system and control
%        theory", https://web.stanford.edu/~boyd/lmibook/lmibook.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       15-March-2021
% Last update:   25-May-2022 (VG, more options, 1D special case)
%                05-July-2022 (VG, removed unecessary input argument)
%                17-March-2023 (VG, bugfix for equally small ellipsoids)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(E(1));

% remove all ellipsoids that only contain a point and add them
q = zeros(n,1);
ind_c = false(size(E));
for i=1:length(E)
    if representsa_(-E(i).q+E(i),'origin',eps)
        q = q +E(i).q;
        ind_c(i) = true;
    end
end
% if only centers, handle
if all(ind_c)
    E = ellipsoid(zeros(n),q);
    return;
end

% at least one "real" ellipsoid remains
E(ind_c) = [];
E(1).q = E(1).q +q;

N = length(E);

% generally, all ellipsoids can be degenerate as long as a positively
% weighted sum of their shape matrices is invertible (otherwise all E_c are
% degenerate and and(...,'outer') will not work)
% Test: If sum of Qi's is positive definite, any sum with non-zero
% weighting is (since Qi's are psd)
Q_sum = zeros(n);
TOL = min(arrayfun(@(ii)E(ii).TOL,length(E)));
for i=1:N
    Q_sum = Q_sum + E(i).Q;
end

T = eye(n);
[V,S,~] = svd(Q_sum);
if all(withinTol(diag(S),0,TOL))
    % result is 0
    E = ellipsoid(zeros(n));
    return;
end
% find max eigenvalue and scaling factor such that max(eig(Q_sum))\approx 1
s = S(1,1);
S = 1/s*S;
for j=1:length(E)
    E(j).Q = 1/s*E(j).Q;
end
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
    % transform all ellipsoids
    E = T'*E;
    for i=1:N
        xt_rem = xt_rem + E(i).q(nt+1:end);
        E(i) = project(E(i),1:nt);
    end
    if n_d==n
        % all Qi's zero
        E = T*ellipsoid(zeros(n),xt_rem);
        return;
    end
end    

% special case: nt==1
% use interval arithmetic to compute exact solution
if nt==1
    Ires = interval(E(1));
    for i=2:length(E)
        Ires = Ires + interval(E(i));
    end
    Et = ellipsoid(rad(Ires)^2,center(Ires));
else
    % if user specified directions, compute sum as intersection of "direction
    % ellipsoids"
    if ~isempty(L)
        % make sure L is normalized
        L = 1./sqrt(sum(L.^2,1)).*L;
    
        E_L = lplus(E,L,mode);
        if strcmp(mode,'outer')
            Et = and_(E_L(1),E_L(2:end),mode);
        elseif strcmp(mode,'inner')
            Et = or(E_L(1),E_L(2:end),mode);
        else
            throw(CORAerror('CORA:wrongValue','third',"'inner' or 'outer'"));
        end
    else
        % choose between fixed point iteration [1] or exact computation 
        % (SDP problem) [2]
        if strcmp(mode,'outer')
            idx = true(length(E),1);
            q = zeros(dim(E(1)),1);
            for j=1:length(E)
                idx(j) = ~representsa_(E(j)-E(j).q,'origin',1e-8);
                q = q + ~idx(j)*E(j).q;
            end
            E = E(idx);
            Et = plusEllipsoidOA(E) + q;
        elseif strcmp(mode,'outer:halder')
            Et = plusEllipsoidOA_halder(E);
        elseif strcmp(mode,'inner')
            throw(CORAerror('CORA:noExactAlg'));
        else
            throw(CORAerror('CORA:wrongValue','third',"'outer' or 'outer:halder'"));
        end
    end
end

% backtransform
E = T'*ellipsoid(s*blkdiag(Et.Q,zeros(n_d)),[Et.q;xt_rem]);

% ------------------------------ END OF CODE ------------------------------
