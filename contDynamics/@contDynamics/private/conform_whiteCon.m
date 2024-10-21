function con = conform_whiteCon(sys,testSuite,p_GO,G_X0,G_U,options_cs,union_y_a)
% conform_whiteCon - formulate reachset conformance constraints as linear
%    constraints
%
% Syntax:
%    con = conform_whiteCon(sys,testSuite,p_GO,G_X0,G_U,options_cs,union_y_a)
%
% Inputs:
%    sys - contDynamics object
%    y_m - cell array containing measurement trajectories of dimension 
%          n_y x N_k for the testcase m
%    p_GO - cell array for mdifferent reference trajectores with the model
%           parameters summarized in a struct with
%               .A{k}: matrix which describes the influence of the initial
%                      state x(1) on the state x(k+1)
%               .B{k,j}: matrix which describes the influence of the input
%                        u(j) on x(k+1)
%               .F{k,j}: matrix which describes the influence of the
%                        linearization error L(j) on x(k+1)  
%               .D{k,j}: matrix which describes the influence of the input
%                        u(j) on the output y(k)
%               .E{k,j}: matrix which describes the influence of the
%                        linearization error L(j) on y(k) 
%               .C{k}: matrix which describes the influence of the initial
%                      state x(1) on y(k)
%               .x: reference state trajectory
%               .u: reference input trajectory
%               .y: reference output trajectory  
%     G_X0 - generators of the initial state set
%     G_U - generators of the disturbance set    
%     options_cs - conformance synthesis specifications
%     union_y_a - n_m x 1 cell array with n_k x n_y x n_s output arrays
%                 with union_y_a{m} = testSuite{m}.y - y_nom
%           
% Outputs:
%     con - struct with equality and inequality constraints, with fields
%              .A_eq: matrix for linear equality contraint
%              .b_eq: vector for linear equality contraint
%              .A_ineq: matrix for linear inequality contraint
%              .b_ineq: vector for linear inequality contraint
%
% References:
%    [1] L. Luetzow and M. Althoff, "Scalable Reachset-conformant 
%        Identification of Linear Systems,"  in IEEE Control Systems
%        Letters, vol. 8, pp. 520-525, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conform_white

% Authors:       Laura Luetzow
% Written:       26-July-2023
% Last update:   16-April-2024 (LL, restructuring)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Q_ineq = [];

if strcmp(options_cs.constraints,"half")
    A_eq = []; b_eq = [];
    [A_ineq,b_ineq] = aux_halfspaceConstraints(sys,testSuite,p_GO,G_X0,G_U,union_y_a);
elseif strcmp(options_cs.constraints,"gen")
    [A_eq,b_eq,A_ineq,b_ineq] = aux_generatorConstraints(testSuite,p_GO,G_X0,G_U,options_cs.n_c);
end

con.A_eq = A_eq;
con.b_eq = b_eq;
con.Q_ineq = Q_ineq;
con.A_ineq = A_ineq;
con.b_ineq = b_ineq;

end


% Auxiliary functions -----------------------------------------------------

function [A_ineq,b_ineq] = aux_halfspaceConstraints(sys, testSuite, ...
    p_GO, G_X0, G_U, union_y_a)

A_ineq = []; b_ineq = [];

for m = 1:length(testSuite)      
    % y_m = testSuite{m}.y;
    n_k = size(p_GO{m}.y, 2); % number of timesteps

    % compute halfspace contraints for each time step
    for k=1:n_k

        % Compute the generators of the reachable output set
        G_Yk = p_GO{m}.C{k}*G_X0;
        for j = 1:k
            G_Yk = [G_Yk p_GO{m}.D{k,j}*G_U];
        end

        % Compute halfspace normal matrix N_k
        Nk = aux_normalVectors(G_Yk);
        if size(Nk,2) == 0
            Nk = zeros(2*size(G_Yk,1),size(G_Yk,1));
            Nk_c = [eye(size(G_Yk,1)); -eye(size(G_Yk,1))];
        else
            Nk_c = Nk;
        end

        % Compute the linear constraint vector q
        sum_D = 0;
        sum_D_Gu_abs = 0;
        for j=1:k
            % Compute sum_j Dj
            sum_D = sum_D + p_GO{m}.D{k,j};

            % Compute sum_j |Nk Dj G_u|
            sum_D_Gu_abs = sum_D_Gu_abs + abs(Nk*p_GO{m}.D{k,j}*G_U);
        end
        % Distance matrix from the uncertain initial state and inputs

        if isa(sys, 'linearSysDT') || isa(sys, 'linearARX')
            A_ineqk = [abs(Nk*p_GO{m}.C{k}*G_X0), sum_D_Gu_abs, Nk_c*p_GO{m}.C{k}, Nk*sum_D];
        else
            A_ineqk = [abs(Nk*p_GO{m}.C{k}*G_X0), sum_D_Gu_abs];
        end
        A_ineq(end+1:end+size(Nk,1),:) = -A_ineqk;

        % compute the linear summand d
        union_y_ak = permute(union_y_a{m}(k,:,:),[2 3 1]);
        b_ineqk = max(Nk_c*union_y_ak,[],2);
        b_ineq(end+1:end+size(Nk,1),1) = -b_ineqk;
    end
end

end

function [A_eq,b_eq,A_ineq,b_ineq] = aux_generatorConstraints(testSuite,p_GO,G_X0,G_U,n_c)

    b_eq = [];

    eta_U = size(G_U,2);
    eta_X0 = size(G_X0,2);
    n_a = eta_X0 + eta_U; % number of scaling factors
    n_y = size(p_GO{1}.C{1},1);
    n_x = size(p_GO{1}.C{1},2);
    n_u = size(p_GO{1}.u,1);
    I_aX0 = [eye(eta_X0), zeros(eta_X0,eta_U)];
    I_aU = [zeros(eta_U,eta_X0), eye(eta_U)];

    Q_b = sparse([]);
    Q_c = sparse([]);
    R_a = sparse([]);

    for m = 1:length(testSuite)
        y_m = testSuite{m}.y;
        n_k = size(p_GO{m}.y, 2); % number of timesteps

        R_as = zeros(n_k * eta_X0 + sum(1:n_k)*eta_U, eta_X0 + eta_U);
        Q_bs = [];
        Q_cs = zeros(n_k*n_y, n_x+n_u);
        I_ak = ([I_aX0; zeros(n_k*eta_U, eta_X0+eta_U)]);
        for k=1:n_k

            % matrices for equality constraints
            % Compute the generators of the reachable output set
            G_Yk = zeros(n_y,eta_X0+k*eta_U);
            G_Yk(:,1:eta_X0) = p_GO{m}.C{k}*G_X0;
            sum_D = 0;
            for j = 1:k
                % Compute sum_j Dj
                if n_c > 0
                    sum_D = sum_D + p_GO{m}.D{k,j};
                end

                G_Yk(:,eta_X0+(j-1)*eta_U+1:eta_X0+j*eta_U) = p_GO{m}.D{k,j}*G_U;
            end
            Q_bs = blkdiag(Q_bs, G_Yk);
            if n_c > 0
                Q_cs((k-1)*n_y+1:k*n_y,:) = [p_GO{m}.C{k} sum_D];
            end

            % matrices for inequality constraints
            I_ak(eta_X0+(k-1)*eta_U+1:eta_X0+k*eta_U,:) = I_aU;
            R_as((k-1)*eta_X0 + sum(1:k-1)*eta_U+1:k*eta_X0 + sum(1:k)*eta_U,:) = I_ak(1:eta_X0+k*eta_U,:);
        end

        % equality constraints Ax = b
        if iscell(y_m)
            n_s = length(y_m); % number of tests
            for s = 1:n_s
                Q_b = blkdiag(Q_b, sparse(Q_bs)); % copy for each sample
                b_eq = [b_eq; reshape(p_GO{m}.y-y_m{s}',[],1)];
            end
        else
            n_s = size(y_m,3);
            Q_b = blkdiag(Q_b, sparse(kron(eye(n_s),Q_bs))); % copy for each sample
            b_eq = [b_eq; reshape(p_GO{m}.y-permute(y_m,[2,1,3]),[],1)];
        end
        if n_c > 0
            Q_c = [Q_c; sparse(-repmat(Q_cs, n_s, 1))];
        end

        % inequality constraints Ax <= b
        R_a = [R_a; sparse(-repmat(R_as, n_s, 1))];
    end
    R_c = sparse(size(R_a,1), n_c);
    R_b = speye(size(R_a,1));

    A_ineq = [R_a R_c R_b; R_a R_c -R_b];
    b_ineq = sparse(size(A_ineq, 1),1);

    A_eq = [sparse(size(Q_b,1), n_a) Q_c Q_b];

    % remove rows which correspond to NaN values due to missing
    % measurements in some test cases
    i_nan = isnan(b_eq); 
    b_eq(i_nan) = [];
    A_eq(i_nan,:) = [];

end

function N = aux_normalVectors(G)
% compute normal vectors of a zonotope with generator matrix G

% remove zero generators
G = nonzeroFilter(G);
if size(G,2) == 0
    N = [];
    return
end
% system dimension and number of generators
[n,nrGen] = size(G);
% get number of possible facets
comb = combinator(nrGen,n-1,'c');
% bypass bug in combinator (rows with all-zeros?!)
comb = comb(any(comb,2),:);

% build matrix of normal vectors
N = zeros(length(comb(:,1)),n);
for i=1:length(comb(:,1))
    indices=comb(i,:);
    Q=G(:,indices);
    v=ndimCross(Q);
    if norm(v) == 0
        N(i,:)=v';
    else
        N(i,:)=v'/norm(v);
    end
end
% second half of normal vectors is symmetric
N = [N; -N];

end

% ------------------------------ END OF CODE ------------------------------
