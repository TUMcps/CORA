function obj = project(obj,dims)
% project - projects a mptPolytope onto a set of dimensions
%
% Syntax:  
%    res = project(obj,dims)
%
% Inputs:
%    obj - mptPolytope object
%    dims - vector of dimensions
%
% Outputs:
%    res - projected mptPolytope object
%
% Example: 
%    poly = mptPolytope.generateRandom(3);
%    poly_ = project(poly,[1,2]);
%
%    figure; hold on; grid on;
%    plot(poly,[1,2,3],'r','FaceAlpha',0.5,'Filled',true);
%    plot(poly_,[1,2],'b','LineWidth',2);
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

    % remove redundant halfspaces
    obj = removeRedundancies(obj);
    
    % determine dimensions that have to be projected away
    n = dim(obj);
    remDims = setdiff(1:n,dims);
    
    % project away all unwanted dimensions by repeatadly applying the
    % Fourier-Motzkin elimination
    for i = 1:length(remDims)
       
        % project away current dimension 
        [A,b] = fourierMotzkinElimination(obj.P.A,obj.P.b,remDims(i));
        
        % update indices to match projected polytope
        remDims = remDims - 1;
        
        % remove redundant halfspaces
        obj = removeRedundancies(mptPolytope(A,b));
    end
    
    % shuffle dimensions of the remaining projected polytope
    [~,ind] = sort(dims);
    A = obj.P.A;
    A(:,ind) = A;
    
    obj = mptPolytope(A,obj.P.b);
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