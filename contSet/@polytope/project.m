function P_out = project(P,dims)
% project - projects a polytope onto a set of dimensions
%
% Syntax:
%    P_out = project(P,dims)
%
% Inputs:
%    P - polytope object
%    dims - vector of dimensions
%
% Outputs:
%    P_out - projected polytope object
%
% Example: 
%    P = polytope.generateRandom('Dimension',3);
%    P_ = project(P,[1,2]);
%
%    figure; hold on;
%    plot(P,[1,2,3],'r','FaceAlpha',0.5);
%    plot(P_,[1,2],'b','LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Authors:       Matthias Althoff, Niklas Kochdumper, Viktor Kotsev
% Written:       14-November-2016
% Last update:   21-December-2020 (NK, use Fourier-Motzkin elimination)
%                28-June-2022 (VK, Changed using the fourier toolbox)
%                15-November-2023 (MW, bug fix for equality constraints)
%                14-July-2024 (MW, support vertex representation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(P);

% no projection, copy the polyhedron
if isequal(dims,1:n)
    P_out = polytope(P);
    return;
end

if any(dims > n)
    throw(CORAerror('CORA:wrongValue','second',...
        sprintf('Cannot compute projection on higher dimension than %i.',n)));
end

% if vertex representation given, much simpler projection
if P.isVRep.val
    P_out = polytope(P.V_.val(dims,:));
    return
end

% default method
method = 'fourier_jones';

% remove zeros
P = compact_(P,'zeros',eps);

% check emptiness
if representsa_(P,'emptySet',eps)
    P_out = polytope.empty(n);
    return;
end

% project the polytope

% remove redundant halfspaces (override input object)
P = compact_(P,'all',1e-9);
% copy input object (rewrite equality constraints as pairwise
% inequality constraints)
P_out = polytope([P.A_.val; P.Ae_.val; -P.Ae_.val],...
    [P.b_.val; P.be_.val; -P.be_.val]);

% determine dimensions that have to be projected away
n = dim(P_out);
removeDims = setdiff(1:n,dims);

% different methods
switch method
    case 'fourier_jones'
        Ab = [P_out.A_.val(:, [dims, removeDims]), P_out.b_.val];
        Ab_ = fourier(Ab, 1:length(dims), 1e-6, 0);
        P_out = polytope(Ab_(:,1:end-1), Ab_(:,end));
        return

    case 'fourier'
        for i = 1:length(removeDims)
           
            % project away current dimension 
            [A,b] = aux_fourierMotzkinElimination(P_out.A_.val,P_out.b_.val,removeDims(i));
            
            % update indices to match projected polytope
            removeDims = removeDims - 1;
            
            % remove redundant halfspaces
            P_out = compact_(polytope(A,b),'all',1e-12);
        end
        
        % sort dimensions of the remaining projected polytope according to dims
        [~,ind] = sort(dims);
        A = P_out.A_.val;
        A(:,ind) = A;
        P_out.A_.val = A;
        return
end

end


% Auxiliary functions -----------------------------------------------------

function [A,b] = aux_fourierMotzkinElimination(A,b,j)
% project the polytope A*x <= b onto the dimension "j" using
% Fourier-Motzkin elimination: algorithm taken from Wikipedia
% https://de.wikipedia.org/wiki/Fourier-Motzkin-Elimination#Die_Fourier-Motzkin-Elimination_aus_Sicht_der_linearen_Algebra

    % number of constraints
    nrCon = size(A,1);

    % divide dim-th column of matrix A into entries = 0, > 0, and < 0
    Z = find(A(:,j)' == 0);
    N = find(A(:,j)' < 0);
    P = find(A(:,j)' > 0);
    
    % compute cartesian product of sets P and N to get all combinations
    list = [kron(N,ones(1,length(P))); repmat(P,[1,length(N)])]';
    nrComb = size(list,1);
    
    % construct projection matrix: number of columns of the projection
    % matrix is the number of constraints of the projected polytope
    m = length(Z);
    U = zeros(m+nrComb,nrCon);

    % for all Z, we have the Z(i)-th basis vector
    for i=1:m
        U(i,Z(i)) = 1;
    end
    
    % for all pairs (s,t) in P x N, we have
    %   a_(t,j) * e_(s) - a_(s,j) * e_(t)
    for i=1:nrComb
        U(m+i,list(i,1)) = A(list(i,2),j);
        U(m+i,list(i,2)) = -A(list(i,1),j);
    end
    % note: resulting matrix U can only have nonnegative entries!

    % perform projection
    A = U * A; A(:,j) = [];
    b = U * b;
end

% ------------------------------ END OF CODE ------------------------------
