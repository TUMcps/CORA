function [A,b,Ae,be] = computeHRep(V)
% computeHRep - conversion from vertex representation to halfspace
%    representation (containing potentially redundant halfspaces)
% 
%    requirements:  the set is full-dimensional
%                   the set has at least d+1 vertices
%
% Syntax:
%    [A,b,Ae,be] = polytope.computeHRep(V)
%
% Inputs:
%    V - Matrix storing vertices ([n,p] for p vertices in dimension n)
%
% Outputs:
%    A,b - inequality constraints: A*x <= b
%    Ae,be - equality constraints: Ae*x = be
%
% Example:
%    V = [1 -1 0 0; 0 0 1 -1];
%    [A,b,Ae,be] = polytope.computeHRep(V)
%
% Reference: MPT-Toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope

% Authors:       Viktor Kotsev
% Written:       06-June-2022
% Last update:   ---
% Last revision: 25-July-2023 (MW, move to private, new return values)

% ------------------------------ BEGIN CODE -------------------------------

% read out dimension and number of vertices
[n,nrVert] = size(V);

% special cases
if nrVert == 1
    % simple method for a single vertex: { x | I * x = V }
    A = zeros(0,n); b = [];
    Ae = eye(n); be = V;
    return

elseif n == 1
    % 1D case
    [A,b,Ae,be] = aux_1D(V);
    return

elseif n == 2 && (nrVert <= 2 || rank(V) < 2)
    % 2D-degenerate case
    [A,b,Ae,be] = aux_2D_degenerate(V);
    return

end

% shift vertices by center
c = mean(V,2); %0.5*(max(V,[],2) + min(V,[],2));
V = V - c;

% check for degeneracy
[Umat,Smat,Vmat] = svd(V);
% dimension of subspace
r = nnz(~withinTol(Smat,0,1e-10));

% project vertices onto subspace
if r < n
    % project vertices onto basis
    V = Umat' * V;
    % filter out subspace
    V = V(1:r,:);
end

% set is actually 1D
if size(V,1) == 1
    A = [1; -1]; 
    b = [max(V); -min(V)];

else
    % compute convex hull (possibly in subspace)
    try
        K = convhulln(V');
    catch ME
        if strcmp(ME.identifier,'MATLAB:qhullmx:UndefinedError')
            % error with precision... take another method
            K = convhulln(V',{'QJ'});
        else
            rethrow(ME);
        end
    end
    
    %interior point
    s.A = zeros(size(K, 1), r);
    s.B = ones(size(K, 1), 1);
    s.lin = [];
    
    for i = 1:size(K, 1)
        % each row of K contains indices of vertices that lie on the i-th
        % facet
        P = V(:, K(i, :))';
        % compute the normal vector and the offset of the facet
        W = [P, -ones(r, 1)];
        [AB, ~] = qr(W'); % qr() is much faster than null()
        a = AB(1:r, end);
        b = AB(end, end);
        
        % determine the sign
        if 0 > b
            a = -a;
            b = -b;
        end
        s.A(i, :) = a';
        s.B(i) = b;
    end
    
    Hall = [s.A s.B];
    
    % inequality representation
    H = Hall; H(s.lin, :) = [];
    A = H(:,1:end-1);
    if r == n
        b = H(:,end) + A*c;
    else
        b = H(:,end);
    end
    
    % equality representation
    He = Hall(s.lin, :);
    Ae = He(:,1:end-1);
    if r == n
        be = He(:,end) + Ae*c;
    else
        be = He(:,end);
    end
end

% back-projection to higher-dimensional space (if computed in subspace)
if r < n
    % extend by zeros
    nrIneq = size(A,1);
    A = [A, zeros(nrIneq,n-r)];
    % back-projection
    A = A * Umat';
    % incorporate original mean into offset
    b = b + A*c;
    
    % equality constraints where singular values are (almost) zero
    Ae = Umat(:,r+1:n)';
    % translate offset by original center
    be = Ae*c;
end

end


% Auxiliary functions -----------------------------------------------------

% function [A,b,Ae,be] = aux_nD_degenerate(V)
% % not enough vertices for standard method -> project into subspace and then
% % project back into higher-dimensional space
% 
% % dimension
% n = size(V,1);
% 
% % shift vertices by center
% c = 0.5*(max(V,[],2) + min(V,[],2)); %mean(V,2);
% V = V - c;
% 
% % find basis
% [Umat,Smat,Vmat] = svd(V);
% % dimension of subspace
% r = nnz(~withinTol(Smat,0,1e-10));
% 
% % project vertices onto basis
% V_ = Umat' * V;
% % filter out subspace
% V_ = V_(1:r,:);
% % compute H-representation in subspace
% [A_,b_,Ae_,be_] = computeHRep(V_);
% 
% % extend by zeros
% nrIneq = size(A_,1);
% A_ = [A_, zeros(nrIneq,n-r)];
% % back-projection
% A = A_ * Umat';
% % incorporate original mean into offset
% b = b_ + A*c;
% 
% % equality constraints where singular values are (almost) zero
% Ae = Umat(:,r+1:n)';
% % translate offset by original center
% be = Ae*c;
% 
% end

function [A,b,Ae,be] = aux_2D_degenerate(V)
% input: V ... vertices ([n,p], n ... dimension, p ... number of points)
% special function to deal with special cases that occur quite often due to
% plotting and other 2D-related functions

Ae = []; be = [];

% only one distinct point given
if size(V,2) == 1 || all(all(withinTol(V - V(:,1),0)))
    % use axis-aligned normal vectors for simplicity (one could also use
    % only three halfspaces, but then other functions become more costly)
    A = [1 0; 0 1; -1 0; 0 -1];
    b = [V(1,1); V(2,1); -V(1,1); -V(2,1)];
    return
end

% collinear vertices -> polytope is just a line

% find minimum and maximum of line using x-value
[~,minIdx] = min(V(1,:));
[~,maxIdx] = max(V(1,:));
% if indices are equal, that means that x-value is constant -> choose y
% value for sorting
if minIdx == maxIdx
    [~,minIdx] = min(V(2,:));
    [~,maxIdx] = max(V(2,:));
end

% compute vector along the line and normalize
dir = V(:,maxIdx) - V(:,minIdx);
dir = dir ./ vecnorm(dir);

% fill in normal vectors, including 90° turn
dir_ = [-dir(2); dir(1)];
A = [dir'; -dir'; dir_'; -dir_'];
b = zeros(4,1);

% compute offset of start and end
b(1) = dir' * V(:,maxIdx);
b(2) = -dir' * V(:,minIdx);
% for orthogonal normal vectors, we can use either maxIdx or minIdx
b(3) = dir_' * V(:,maxIdx);
b(4) = -dir_' * V(:,maxIdx);

end

function [A,b,Ae,be] = aux_1D(V)
% input: V ... vertices ([n,p], n ... dimension, p ... number of points)

% find minimum and maximum value
maxV = max(V);
minV = min(V);
% check unboundedness
maxV_inf = maxV == Inf;
minV_inf = minV == -Inf;

% generally two inequality constraints
% x <= maximum value && x >= minimum value
if maxV_inf && minV_inf
    % both are infinity -> no constraints
    A = []; b = [];
    Ae = []; be = [];

elseif maxV_inf && ~minV_inf
    % only bounded from below
    A = -1; b = -minV;
    Ae = []; be = [];

elseif ~maxV_inf && minV_inf
    % only bounded from above
    A = 1; b = maxV;
    Ae = []; be = [];

elseif withinTol(maxV,minV)
    % single point -> use equality constraint
    A = []; b = [];
    Ae = 1; be = maxV;

else
    % bounded from above and below, not a single point
    A = [1; -1]; b = [maxV; -minV];
    Ae = []; be = [];
end

end

% ------------------------------ END OF CODE ------------------------------
