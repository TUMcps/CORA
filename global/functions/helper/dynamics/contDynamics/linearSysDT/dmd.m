function [A,B,c] = dmd(X,dt,varargin)
% dmd - perform Dynamic Mode Decomposition (DMD) to approximate the data
%    points by a discrete-time linear system
%
% Syntax:
%    [A,B,c] = dmd(X,dt)
%    [A,B,c] = dmd(X,dt,U)
%    [A,B,c] = dmd(X,dt,U,type)
%    [A,B,c] = dmd(X,dt,U,type,alg)
%    [A,B,c] = dmd(X,dt,U,type,'svd',rank)
%
% Inputs:
%    X - data matrix, where the number of columns is equal to the number of
%        data points. For multiple traces X can be specified as cell-array
%    dt - time step size of the data
%    U - inputs corresponding to the data matrix, where the number of
%        colums is equal to the number of data points minus one 
%    type - return discrete time system (type = 'disc') or continuous time
%           system (type = 'cont'). The default value is type = 'disc'.
%    alg - algorithm used for DMD ('svd' or 'arnoldi'). The default value 
%          is 'svd'
%    rank - rank of the truncated singular value decomposition
%
% Outputs:
%    A - system matrix for the appox. discrete time system x[i+1] = A*x[i]
%    B - input marix for the approx. system x[i+1] = A*x[i] + B*u[i]
%    c - constant input for the approx. system x[i+1] = A*x[i] + B*u[i] + c
%
% Example: 
%    % discrete-time linear system
%    A = [-0.4 0.6; 0.6 -0.4]; B = [0;1]; dt = 0.1;
%    sys = linearSysDT(A,B,dt);
%
%    % generate data points by simulating the disturbed system
%    params.R0 = interval([0.9;1.9],[1;2]);
%    params.U = interval(-0.1,0.1);
%    params.tFinal = 1;
%    options.points = 1;
%    simRes = simulateRandom(sys,params,options);
%    X = simRes.x{1}';
%
%    % use dynamic mode decomposition to estimate the matrix A from data
%    A_ = dmd(X)
%    A
%
% References:
%    [1] J. Proctor and et al. "Including inputs and control within
%        equation-free architectures for complex systems", 2016
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       17-November-2021 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % default values
    [U,type,alg,rank] = setDefaultValues({[],'disc','svd',[]},varargin);

    % check input arguments
    inputArgsCheck({{type,'str',{'disc','cont'}},...
        {alg,'str',{'svd','arnoldi'}}});

    % divide data into original matrix X1 and time-shifted matrix X2
    if iscell(X)
        X1 = []; X2 = [];
        for i = 1:length(X)
            if strcmp(type,'cont')
                X1 = [X1 X{i}(:,2:end-1)];
                X2 = [X2 (X{i}(:,3:end) - X{i}(:,1:end-2))./(2*dt)];
            else
                X1 = [X1 X{i}(:,1:end-1)]; X2 = [X2 X{i}(:,2:end)]; 
            end
        end
    else
        if strcmp(type,'cont')
            X1 = X(:,2:end-1); X2 = (X(:,3:end) - X(:,1:end-2))./(2*dt);
        else
            X1 = X(:,1:end-1); X2 = X(:,2:end);
        end
    end
    
    % add inputs (see Sec. 3.1 in [1])
    if ~isempty(U)
        if iscell(U)
           U_ = [];
           for i = 1:length(U)
              if strcmp(type,'cont')
                 U_ = [U_,(U{i}(:,1:end-1)+U{i}(:,2:end))/2];  
              else
                 U_ = [U_,U{i}];
              end
           end
           X1 = [X1;U_];
        else
            if strcmp(type,'cont')
                X1 = [X1; (U(:,1:end-1)+U(:,2:end))/2];
            else
                X1 = [X1;U];
            end
        end
    end
    
    % apply different algorithms
    if strcmp(alg,'svd')

        % singular value decomposition
        [V,S,W] = svd(X1,'econ');
        
        % reduce rank by removing the smallest singular values
        if ~isempty(rank) && rank < size(S,1)
            V = V(:,1:rank); S = S(1:rank,1:rank); W = W(:,1:rank);
        end
        
        % compute resulting system matrix A (see Eq. (20) in [1])
        A = X2*W*diag(1./diag(S))*V';
    
    elseif strcmp(alg,'arnoldi')
        
        % compute coefficient vector a using QR-decomposition
        [Q,R] = qr(X2);
        a = linsolve(R,Q'*X(:,end));
        
        % perform an eigenvalue decomposition of the companion matrix S
        S = diag(ones(length(a)-1,1),-1);
        S(:,end) = a;
        
        [V,D] = eig(S);
        
        % select the largest eigenvalues
        [~,ind] = sort(abs(diag(D)),'descend');
        ind = ind(1:size(X,1));
        
        D = D(ind,ind); V = V(:,ind);
        
        % compute approximate eigenvectors
        V = X1*V;
        
        % construct system matrix A
        A = real(V*D*inv(V));  
    end
    
    % devide into state matrix and input matrix
    if ~isempty(U)
        B = A(:,size(A,1)+1:end); A = A(:,1:size(A,1));
    else
        B = [];
    end
    
    % determine constant offset c
    if nargout > 2
        c = aux_compConstOffset(A,B,X,U,dt,type);
    end
end


% Auxiliary functions -----------------------------------------------------

function c = aux_compConstOffset(A,B,X,U,dt,type)
% compute the optimal constant offset c to match the data
    
    % create all-zero inputs if inputs are empty
    if isempty(B)
       B = zeros(size(A,1),1);
       U = cell(length(X),1);
       for i = 1:length(U)
           U{i} = zeros(1,size(X{i},2)-1); 
       end
    end

    % get propagation matrices
    F = eye(size(A,1));
    if strcmp(type,'cont')
        sys = linearSysDT(linearSys(A,F),dt);
        A = sys.A; B = sys.B*B; F = sys.B;
    end

    % construct data matrices X1*c = X2
    X1 = []; X2 = [];
    
    for i = 1:length(X)

       X1_ = []; X2_ = [];
       tmp = F; pred = A*X{i}(:,1) + B*U{i}(:,1); 

       for j = 2:size(X{i},2)
          X1_ = [X1_; tmp]; tmp = A*tmp + F;
          X2_ = [X2_; X{i}(:,j) - pred]; 
          if j < size(X{i},2)
              pred = A*pred + B*U{i}(:,j);
          end
       end
    end

    % estimate pseudo-inverse using a singular value decomposition
    [U,S,W] = svd(X1,'econ');
    c = W*diag(1./diag(S))*U'*X2;
end
    
% ------------------------------ END OF CODE ------------------------------
