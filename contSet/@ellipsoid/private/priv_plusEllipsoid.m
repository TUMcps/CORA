function E = priv_plusEllipsoid(E_cell,L,mode)
% priv_plusEllipsoid - Computes an inner approximation or outer
%    approximation of the Minkowski sum of a list of ellipsoids
%
% Syntax:
%    E = priv_plusEllipsoid(E_cell,L,mode)
%
% Inputs:
%    E_cell - cell array of ellipsoid objects
%    L - directions 
%    mode - type of approximation: 'inner', 'outer'
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
% See also: ellipsoid/plus

% Authors:       Victor Gassmann
% Written:       15-March-2021
% Last update:   25-May-2022 (VG, more options, 1D special case)
%                05-July-2022 (VG, removed unecessary input argument)
%                17-March-2023 (VG, bugfix for equally small ellipsoids)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension
n = dim(E_cell{1});

% check which ellipsoids only contain a single point
idx_isPoint = cellfun(@(E_i) representsa_(E_i,'point',eps),E_cell,'UniformOutput',true);
q = zeros(n,1);
if any(idx_isPoint)
    q = sum(cell2mat(cellfun(@(E_i) E_i.q,E_cell(idx_isPoint),'UniformOutput',false)),2);
end

% if all consist only of the center, exit
if all(idx_isPoint)
    E = ellipsoid(zeros(n),q);
    return;
end

% remove all ellipsoid that only contain a point
E_cell(idx_isPoint) = [];
% add the accumulated shift to the first remaining ellipsoid
E_cell{1}.q = E_cell{1}.q + q;
% count the number of remaining ellipsoids
N = length(E_cell);

% if only one remains, we can also exit
if N == 1
    E = E_cell{1};
    return
end

% generally, all ellipsoids can be degenerate as long as a positively
% weighted sum of their shape matrices is invertible (otherwise all E_c are
% degenerate and and(...,'outer') will not work)
% Test: If sum of Qi's is positive definite, any sum with non-zero
% weighting is (since Qi's are psd)
Q_sum = sum(reshape(cell2mat(cellfun(@(E_i) E_i.Q,E_cell,'UniformOutput',false)),n,n,N),3);
TOL = min(cellfun(@(E_i) E_i.TOL,E_cell,'UniformOutput',true));

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
for j=1:length(E_cell)
    E_cell{j}.Q = 1/s*E_cell{j}.Q;
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
    E_cell = T'*E_cell;
    for i=1:N
        xt_rem = xt_rem + E_cell{i}.q(nt+1:end);
        E_cell{i} = project(E_cell{i},1:nt);
    end
    if n_d == n
        % all Qi's zero
        E = T*ellipsoid(zeros(n),xt_rem);
        return;
    end
end    

% special case: nt==1
% use interval arithmetic to compute exact solution
if nt == 1
    Ires = interval(E_cell{1});
    for i=2:length(E_cell)
        Ires = Ires + interval(E_cell{i});
    end
    Et = ellipsoid(rad(Ires)^2,center(Ires));
else
    % if user specified directions, compute sum as intersection of "direction
    % ellipsoids"
    if ~isempty(L)
        % make sure L is normalized
        L = 1./sqrt(sum(L.^2,1)).*L;
    
        E_L = priv_lplus(E_cell,L,mode);
        if strcmp(mode,'outer')
            Et = and_(E_L(1),E_L(2:end),mode);
        elseif strcmp(mode,'inner')
            Et = or(E_L(1),E_L(2:end),mode);
        else
            throw(CORAerror('CORA:wrongValue','third',...
                "'inner' or 'outer'"));
        end
    else
        % choose between fixed point iteration [1] or exact computation 
        % (SDP problem) [2]
        if strcmp(mode,'outer')
            idx = true(length(E_cell),1);
            q = zeros(dim(E_cell{1}),1);
            for j=1:length(E_cell)
                idx(j) = ~representsa_(E_cell{j},'point',1e-8);
                q = q + ~idx(j)*E_cell{j}.q;
            end
            E_cell = E_cell(idx);
            Et = priv_plusEllipsoidOA(E_cell) + q;
        elseif strcmp(mode,'outer:halder')
            Et = priv_plusEllipsoidOA_halder(E_cell);
        elseif strcmp(mode,'inner')
            throw(CORAerror('CORA:noExactAlg'));
        else
            throw(CORAerror('CORA:wrongValue','third',...
                "'outer' or 'outer:halder'"));
        end
    end
end

% backtransform
E = T'*ellipsoid(s*blkdiag(Et.Q,zeros(n_d)),[Et.q;xt_rem]);

% ------------------------------ END OF CODE ------------------------------
