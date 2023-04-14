function P = project(P,dims)
% project - projects a polytope onto the specified dimensions
%
% Syntax:  
%    res = project(P,dims)
%
% Inputs:
%    P - (mptPolytope) polytope
%    dims - dimensions for projection
%
% Outputs:
%    res - (mptPolytope) projected polytope
%
% Example: 
%    P = mptPolytope([1 0 0; 0 1 0; 0 0 1; -1 -1 -1],[1;1;1;1]);
%    P_ = project(P,[1,2]);
%
%    figure; hold on; grid on;
%    plot(P,[1,2,3],'FaceColor','r','FaceAlpha',0.5);
%    plot(P_,[1,2],'b','LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-November-2016
% Last update:  21-December-2020 (NK, use Fourier-Motzkin elimination)
% Last revision:---

%------------- BEGIN CODE --------------

    % if dimensions given as logical array, convert to numbers
    if islogical(dims)
        temp = 1:length(dims);
        dims = temp(dims);
    end

    % remove redundant halfspaces
    P = removeRedundancies(P);
    
    % determine dimensions that have to be projected away
    n = dim(P);
    remDims = setdiff(1:n,dims);
    
    % project away all unwanted dimensions by repeatadly applying the
    % Fourier-Motzkin elimination
    for i = 1:length(remDims)
       
        % project away current dimension 
        [A,b] = fourierMotzkinElimination(P.P.A,P.P.b,remDims(i));
        
        % update indices to match projected polytope
        remDims = remDims - 1;
        
        % remove redundant halfspaces
        P = removeRedundancies(mptPolytope(A,b));
    end
    
    % shuffle dimensions of the remaining projected polytope
    [~,ind] = sort(dims);
    A = P.P.A;
    A(:,ind) = A;
    
    P = mptPolytope(A,P.P.b);
end


% Auxiliary Function ------------------------------------------------------

function [A,b] = fourierMotzkinElimination(A,b,dims)
% project the polytope A x <= b to the dimension "dim" using
% Fourier-Motzkin elimination

    % divide dim-th column of matrix A into entries = 0, > 0, and < 0
    Z = find(A(:,dims)' == 0);
    N = find(A(:,dims)' < 0);
    P = find(A(:,dims)' > 0);
    
    % compute cartesian product of sets P and N to get all combinations
    list = [kron(N,ones(1,length(P))); repmat(P,[1,length(N)])];
    
    % construct projection matrix
    m = length(Z);
    U = zeros(m + size(list,2),size(A,1));

    for i = 1:length(Z)
       U(i,Z(i)) = 1;
    end
    
    for i = 1:size(list,2)
       U(m+i,list(1,i)) = A(list(2,i),dims);
       U(m+i,list(2,i)) = -A(list(1,i),dims);
    end

    % projection
    A = U * A; A(:,dims) = [];
    b = U * b;
end

%------------- END OF CODE --------------