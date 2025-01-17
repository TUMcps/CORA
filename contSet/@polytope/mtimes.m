function P_out = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
%    a polytope
%
% Syntax:
%    P_out = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix/interval matrix/polytope object
%    factor2 - numerical matrix/interval matrix/polytope object
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
%    M = [2 1; -1 2];
%    P_mtimes = M*P;
% 
% Reference: 
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Viktor Kotsev, Mark Wetzlinger
% Written:       01-February-2011
% Last update:   10-September-2015
%                15-June-2016
%                25-July-2016 (intervalhull replaced by interval)
%                28-June-2022
%                12-September-2022 (add affineMap/invAffineMap functions)
%                14-November-2023 (MW, handling for V-rep, scaling method)
%                04-March-2024 (TL, allow right multiplication with scalar)
%                14-July-2024 (MW, implement projective multiplication)
%                28-August-2024 (MW, implement lifting multiplication)
%                29-November-2024 (MW, multiplication with all-zero matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check dimensions
equalDimCheck(factor1,factor2);

% order arguments correctly
[P_copy,matrix] = findClassArg(factor1,factor2,'polytope');
% copy polytope
P_out = polytope(P_copy);

% read out dimension
n = dim(P_out);

% fullspace
if representsa_(P_out,'fullspace',0)
    P_out = polytope.Inf(n);
    return;
end

% numeric matrix
if isnumeric(matrix)

    % map with all-zero matrix
    if all(withinTol(matrix,0,1e-12),'all')
        P_out = polytope.origin(n);
        return
    end

    if isa(factor1,'polytope')
        if isscalar(matrix)
            P_out = matrix * P_out;
            return
        else
            % polytope * matrix case... definition?
            throw(CORAerror('CORA:noops',factor1,factor2));
        end
    end

    % special method for scaling only (and 1D)
    if length(matrix) == 1
        P_out = aux_mtimes_scaling(P_out,matrix);
        return
    end

    % simple method if matrix is invertible
    if diff(size(matrix)) == 0 && rank(matrix) == length(matrix)
        P_out = aux_mtimes_square_inv(P_out,matrix);
        return
    end

    % quicker computation using V-representation
    if P_out.isVRep.val
        P_out = polytope(matrix*P_out.V);
    else
        % method for general mappings
        [m,n] = size(matrix);
        if m < n
            P_out = aux_mtimes_projection(P_out,matrix);
        else
            P_out = aux_mtimes_lifting(P_out,matrix);
        end
    end
    
elseif isa(matrix,'interval') || isa(matrix,'intervalMatrix')

    % only supported for square matrices
    if diff(size(matrix)) ~= 0
        throw(CORAerror('CORA:notSupported',['Multiplication of interval ' ...
            'matrix with polytope only supported for square interval matrices.']));
    end

    % get minimum and maximum
    M_min = infimum(matrix);
    M_max = supremum(matrix);
    % get center of interval matrix
    T = 0.5*(M_max+M_min);
    % get symmetric interval matrix
    Sval = 0.5*(M_max-M_min);
    S = interval(-Sval,Sval);

    % compute interval of polytope
    I = interval(P_out);

    % polytope of interval computations
    Iadd = S*I;
    Padd = polytope(Iadd);

    % compute new polytope
    P_out = T*P_out + Padd;

else
    % specifically, matZonotope and matPolytope multiplication is not
    % supported
    throw(CORAerror('CORA:noops',factor1,factor2));

end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_mtimes_scaling(P,fac)
% simple method for scaling:
%    M S = { M s | s in S }
% -> fac * I * S = { fac * I * x | A x <= b, Ae x == b}    set y = fac*I*x
%                = { y | A (1/fac*I*y) <= b, Ae (1/fac*I*y) == b }
%                = { y | A/fac y <= b, Ae/fac y == b }
%                = { y | A y <= b*fac, Ae y == b*fac }     if fac > 0
%             OR = { y | -A y <= -b*fac, Ae y == b*fac }   if fac < 0
% (note: case with fac = 0 yields a polytope that is just the origin)

if fac == 0
    % resulting polytope is only the origin
    P = polytope.origin(dim(P));
    return
elseif fac > 0
    P.b_.val = P.b_.val * fac;
else
    P.A_.val = -P.A_.val;
    P.b_.val = -P.b_.val * fac;
end
P.be_.val = P.be_.val * fac;

% map vertices if given
if P.isVRep.val
    P.V_.val = P.V_.val * fac;
end

end

function P = aux_mtimes_square_inv(P,M)
% matrix M is invertible

% apply well-known formula (constraints times inverse of matrix)
if P.isHRep.val
    P.A_.val = P.A_.val / M;
    P.Ae_.val = P.Ae_.val / M;
end

% map vertices if given
if P.isVRep.val
    P.V_.val = M*P.V_.val;
end

end

function P = aux_mtimes_projection(P,M)
% general matrix multiplication with projection, see [1, (24)]

n = dim(P);

% compute singular value decomposition of the matrix
[U,S,V] = svd(M);
% number of singular values
r = nnz(S > 1e-10);
% get inverted diagonal matrix with non-zero singular values
D_inv = diag(1./diag(S(1:r,1:r)));

% init polytope before projection
P.A_.val = P.A_.val * V * [D_inv zeros(r,n-r); zeros(n-r,r) eye(n-r)];
P.Ae_.val = P.Ae_.val * V * [D_inv zeros(r,n-r); zeros(n-r,r) eye(n-r)];

% project onto first r dimensions
P = project(P,1:r);

% multiply with orthogonal matrix
P = U*P;

end

function P = aux_mtimes_lifting(P,M)

% read out dimensions
n = dim(P);
m = size(M,1);

% number of inequality constraints and equality constraints
nrIneq = length(P.b_.val);
nrEq = length(P.be_.val);

% compute singular value decomposition of the matrix
[U,S,V] = svd(M);
% number of singular values
r = nnz(S > 1e-10);
% get inverted diagonal matrix with non-zero singular values
D_inv = diag(1./diag(S(1:r,1:r)));

% init polytope before projection
P.A_.val = [P.A_.val * V * D_inv, zeros(nrIneq,m-r)] * U';
P.Ae_.val = [P.Ae_.val * V * D_inv, zeros(nrEq,m-r); zeros(m-r,n), eye(m-r)] * U';

P.be_.val = [P.be_.val; zeros(m-r,1)];

end

% ------------------------------ END OF CODE ------------------------------
