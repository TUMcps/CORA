function P = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
%    a polytope
%
% Syntax:
%    P = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix/interval matrix/polytope object
%    factor2 - numerical matrix/interval matrix/polytope object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    ---
% 
% Reference: MPT-Toolbox https://www.mpt3.org/
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check dimensions
equalDimCheck(factor1,factor2);

% order arguments correctly
[P,matrix] = findClassArg(factor1,factor2,'polytope');

% empty polytope
if isemptyobject(P)
    P = polytope();
    return;
end

%get dimension
n = dim(P);

%numeric matrix
if isnumeric(matrix)

    % quick computation using V-representation
    if ~isempty(P.V.val)
        if length(matrix)==1
            matrix = matrix*eye(n);
        end
        P = polytope(matrix*P.V.val);
        return
    end

    % special method for scaling
    if length(matrix)==1
        if matrix == 0
            P = polytope(zeros(n,1));
        else
            P = aux_scaling(P,matrix);
        end
    elseif isa(factor1,'polytope')
        P = aux_invAffineMap(P,matrix);
    else
        P = aux_affineMap(P,matrix);
    end
    
else
    %interval matrix 
    if isa(matrix,'interval')
        %get minimum and maximum
        M_min=infimum(matrix);
        M_max=supremum(matrix);
    elseif isa(matrix,'intervalMatrix')
        %get minimum and maximum
        M_min=infimum(matrix.int);
        M_max=supremum(matrix.int);
    end
    %get center of interval matrix
    T=0.5*(M_max+M_min);
    %get symmetric interval matrix
    Sval = 0.5*(M_max-M_min);
    S = interval(-Sval,Sval);

    %compute interval of polytope
    I = interval(P);

    %polytope of interval computations
    Iadd = S*I;
    Padd = polytope(Iadd);

    %compute new polytope
    P=T*P + Padd; 
end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_scaling(P,fac)
% simple method for scaling:
%    M S = { M s | s in S }
% -> fac * I * S = { fac * I * x | A x <= b, Ae x == b}    set y = fac*I*x
%                = { y | A (1/fac*I*y) <= b, Ae (1/fac*I*y) == b }
%                = { y | A/fac y <= b, Ae/fac y == b }
%                = { y | A y <= b*fac, Ae y == b*fac }     if fac > 0
%             OR = { y | -A y <= -b*fac, Ae y == b*fac }   if fac < 0
% (note: case with fac = 0 handled outside)

if fac > 0
    P.b = P.b * fac;
else
    P.A = -P.A;
    P.b = -P.b * fac;
end
P.be = P.be * fac;

end

function P = aux_affineMap(P,M)
% affineMap - computes the affine map M*P of the polytope P with the matrix
%    M of the size n x d.
% 
%    If n  < d then this is projection
%    If n  > d then this is a lifting
%    If n == d then this is rotation/skew.
%
% Syntax:  
%    P = affineMap(P,M)
%
% Inputs:
%    P - polytope object
%    M - matrix
%
% Outputs:
%    P - polytope object

% emptiness checked elsewhere

% if norm(M) <= 1E-12
% 	% Special case: zero map. Return a singleton (the origin) of
% 	% appropriate dimension: M*P = { z | z=0 } where the dimension of "z"
% 	% is equal to the number of rows of "M"
% 	new_dim = size(M, 1);
% 	V = zeros(new_dim, 1);
% 	He = [eye(new_dim), zeros(new_dim, 1)];
% 	P = polytope([], [], He(:,end-1), He(:,end));
% 	return
% end


if ~isempty(P.Ae)
    if isempty(P.A)
        % edge case: no inequalities
        if size(M, 1)==size(M, 2) && abs(det(M)) > 1E-12
            % simple solution if "M" is invertible
            P = polytope([],[],P.Ae*inv(M),P.be);
            return
        else
            % lower-dimensional mapping, solution as suggester by Magnus
            % Nilsson
            nM = size(M,1);
            p = P.Ae\P.be; % particular solution to Ae*d=be
            Ae = zeros(nM); % a bit ugly way to avoid empty matrices for Ae when it should be zero matrices.
            Ae = [null([M*null(P.Ae)]')';Ae]; % Concatenate with the zero matrix...
            Ae = Ae(1:nM, :);                 % ...and remove superfluous rows.
            be = Ae*M*p;
            P = polytope([],[],Ae,be);
            return
        end
        
    elseif size(M,1) == size(M,2)
        % rewrite equalities as pairwise inequalities
        P = polytope([P.A; P.Ae; -P.Ae],[P.b; P.be; -P.be]);

    else
        throw(CORAerror('CORA:notSupported','Case currently not supported.'));

    end
end

% Compute permutation of M s.t. P*y = [M1 M2; M3 M4]*[xr;xn] with rank M1 = rank T and Q*x=[xr;xn]
[L,U,p,q] = lu(sparse(M),'vector');
r = rank(M,1e-12);
pr = p(1:r); pn = p(r+1:end);
qr = q(1:r); qn = q(r+1:end);

rk = rank(full(U(:,1:r)));
if rk~=r
    % if invertibility is not achieved try reduced echelon elimination
    [~,jb]=rref(M,1E-12);
    [L,U,p] = lu(sparse(M(:,jb)),'vector');
    % update column selection
    q = [jb, setdiff(1:dim(P),jb)];
    pr = p(1:r); pn = p(r+1:end);
    qr = q(1:r); qn = q(r+1:end);
end

beta = M(pr,qr) \ eye(r);
 
% Build polytope who's projection is the full-dimensional portion of the mapping
A = P.A;
Ptmp = polytope(A(:,q)*[beta -beta*M(pr,qn); zeros(length(qn),r) eye(length(qn))], P.b);
ptmp = project(Ptmp,1:r);
 
A  = zeros(size(ptmp.A,1),size(M,1));
Ae = zeros(length(pn),size(M,1));
A(:,pr) = ptmp.A;
b = ptmp.b;
Ae(:,pr) = M(pn,qr)*beta;
Ae(:,pn) = -eye(length(pn));
be = zeros(size(Ae,1),1);

% instantiate polytope
P = polytope(A, b, Ae, be);

end

function P = aux_invAffineMap(P,T,t)
% computes P*M
%
% Syntax:  
%    P = invAffineMap(P, T)
%	 P = invAffineMap(P, T, t)
%
% Inputs:
%    P - polytope object  (n x m)
%    T - Square matrix (m x m)
%	 t - vector (n x 1)
%
% Outputs:
%    P - polytope object

if nargin<3
    t = zeros(dim(P(1)), 1);
end

% Only supports square mapping
if size(T, 1)~=size(T, 2)
	error('Only square mappings supported.');
end

if isempty(P.Ae)
	% faster call if we have no equalities
	P = polytope(P.A*T, P.b-P.A*t);
else
	P = polytope(P.A*T, P.b-P.A*t, P.Ae*T, P.be - P.Ae*t);
end

end


% ------------------------------ END OF CODE ------------------------------
