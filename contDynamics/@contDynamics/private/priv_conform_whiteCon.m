function con = priv_conform_whiteCon(sys,testSuite,p_GO,G_X0,G_U,options_cs,union_y_a,n_c)
% priv_conform_whiteCon - formulate reachset conformance constraints as linear
%    constraints
%
% Syntax:
%    con = priv_conform_whiteCon(sys,testSuite,p_GO,G_X0,G_U,options_cs,union_y_a)
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
%     union_y_a - n_m x 1 cell array with dim_y x n_k x n_s output arrays
%                 with union_y_a{m} = testSuite(m).y - y_nom
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
% See also: contDynamics/priv_conform_white

% Authors:       Laura Luetzow
% Written:       26-July-2023
% Last update:   16-April-2024 (LL, restructuring)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(options_cs.constraints,"half")
    A_eq = []; b_eq = [];
    [A_ineq,b_ineq] = aux_halfspaceConstraints(sys,testSuite,p_GO,G_X0,G_U,union_y_a,n_c,options_cs);
elseif strcmp(options_cs.constraints,"gen")
    [A_eq,b_eq,A_ineq,b_ineq] = aux_generatorConstraints(sys,testSuite,p_GO,G_X0,G_U,n_c,options_cs);
end

con.A_eq = A_eq;
con.b_eq = b_eq;
con.A_ineq = A_ineq;
con.b_ineq = b_ineq;

end


% Auxiliary functions -----------------------------------------------------

function [A_ineq,b_ineq] = aux_halfspaceConstraints(sys,testSuite, ...
    p_GO, G_X0, G_U, union_y_a,n_c,options_cs)

A_ineq = []; b_ineq = [];
A_ineq_z = []; A_ineq_T = []; % only required for classification tasks
n_y = size(testSuite(1).y,1);

for m = 1 : length(testSuite)
    % y_m = testSuite(m).y;
    if isa(sys, 'linearARX') || isa(sys, 'linearSysDT')
        % do only one iteration for linear systems
        p_GO_m = p_GO;
        union_y_a_m = union_y_a;
        n_k = size(union_y_a,2);
        if m > 1
            continue
        end
    else
        p_GO_m = p_GO{m};
        union_y_a_m = union_y_a{m};
        n_k = testSuite(m).n_k;
    end

    % compute halfspace contraints for each time step
    for k=1:n_k

        % Compute the generators of the reachable output set
        G_Yk = p_GO_m.C{k}*G_X0;
        for j = 1:k
            G_Yk = [G_Yk p_GO_m.D{k,j}*G_U];
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
            sum_D = sum_D + p_GO_m.D{k,j};

            % Compute sum_j |Nk Dj G_u|
            sum_D_Gu_abs = sum_D_Gu_abs + abs(Nk*p_GO_m.D{k,j}*G_U);
        end
        % Distance matrix from the uncertain initial state and inputs

        if n_c > 0
            A_ineqk = [abs(Nk*p_GO_m.C{k}*G_X0), sum_D_Gu_abs, Nk_c*p_GO_m.C{k}, Nk*sum_D];
        else
            A_ineqk = [abs(Nk*p_GO_m.C{k}*G_X0), sum_D_Gu_abs];
        end
        A_ineq(end+1:end+size(Nk,1),:) = -A_ineqk;

        if strcmp(options_cs.task,"class")
            b_ineqk = Nk_c*p_GO_m.y(:,k);
            b_ineq(end+1:end+size(Nk,1),1) = b_ineqk;

            A_ineq_z = blkdiag(A_ineq_z,Nk_c);
            T_is = [];
            for i=1:n_y
                if testSuite(m).y(i,1) ~= 1 % false class
                    continue
                end
                T_is = [T_is; repmat((1:n_y)==i, n_y, 1) - eye(n_y)];
            end
            if isempty(T_is)
                T_is = zeros(1,size(Nk_c,2));
            end
            A_ineq_T = blkdiag(A_ineq_T, -T_is);
        else
            % compute the linear summand d
            union_y_ak = permute(union_y_a_m(:,k,:),[1 3 2]);
            b_ineqk = max(Nk_c*union_y_ak,[],2);
            b_ineq(end+1:end+size(Nk,1),1) = -b_ineqk;
        end
    end
end

if strcmp(options_cs.task,"class")
    A_ineq = [A_ineq A_ineq_z; ...
        zeros(size(A_ineq_T,1),size(A_ineq,2)) A_ineq_T];
    b_ineq = [b_ineq; zeros(size(A_ineq_T,1),1)];
end
end

function [A_eq,b_eq,A_ineq,b_ineq] = aux_generatorConstraints(sys,testSuite,p_GO,G_X0,G_U,n_c,options_cs)

b_eq = [];

eta_U = size(G_U,2);
eta_X0 = size(G_X0,2);
n_a = eta_X0 + eta_U; % number of scaling factors
I_aX0 = [eye(eta_X0), zeros(eta_X0,eta_U)];
I_aU = [zeros(eta_U,eta_X0), eye(eta_U)];

Q_b = sparse([]);
Q_c = sparse([]);
R_a = sparse([]);

for m = 1:length(testSuite)
    if isa(sys, 'linearARX') || isa(sys, 'linearSysDT')
        p_GO_m = p_GO; % same GO parameters for all trajectories
        y_p_GO_m = p_GO_m.y(:,:,m);
    else
        p_GO_m = p_GO{m};
        y_p_GO_m = p_GO_m.y;
    end
    n_y = size(p_GO_m.C{1},1);
    n_x = size(p_GO_m.C{1},2);
    n_u = size(p_GO_m.u,1);
    y_m = testSuite(m).y;
    n_k = size(y_m, 2); % number of timesteps

    T_i = eye(n_y);
    if strcmp(options_cs.task,"class")
        n_classes = n_y; % loop through all classes
    else
        n_classes = 1;
    end
    for i=1:n_classes
        if strcmp(options_cs.task,"class") % classification
            if y_m(i,1) ~= 1 % false class
                continue
            end
            T_i = repmat((1:n_y)==i, n_y, 1) - eye(n_y);
        end
        R_as = zeros(n_k * eta_X0 + sum(1:n_k)*eta_U, eta_X0 + eta_U);
        Q_bs = [];
        Q_cs = zeros(n_k*n_y, n_x+n_u);
        I_ak = ([I_aX0; zeros(n_k*eta_U, eta_X0+eta_U)]);
        for k=1:n_k

            % matrices for equality constraints
            % Compute the generators of the reachable output set
            G_Yk = zeros(n_y,eta_X0+k*eta_U);
            G_Yk(:,1:eta_X0) = p_GO_m.C{k}*G_X0;
            sum_D = 0;
            for j = 1:k
                % Compute sum_j Dj
                if n_c > 0
                    sum_D = sum_D + p_GO_m.D{k,j};
                end

                G_Yk(:,eta_X0+(j-1)*eta_U+1:eta_X0+j*eta_U) = p_GO_m.D{k,j}*G_U;
            end
            Q_bs = blkdiag(Q_bs, T_i*G_Yk);
            if n_c > 0
                Q_cs((k-1)*n_y+1:k*n_y,:) = [p_GO_m.C{k} sum_D];
            end

            % matrices for inequality constraints
            I_ak(eta_X0+(k-1)*eta_U+1:eta_X0+k*eta_U,:) = I_aU;
            R_as((k-1)*eta_X0 + sum(1:k-1)*eta_U+1:k*eta_X0 + sum(1:k)*eta_U,:) = I_ak(1:eta_X0+k*eta_U,:);
        end

        % equality constraints Ax = b
        n_s = size(y_m,3);
        if strcmp(options_cs.task,"class") % classification
            b_eq = [b_eq; T_i*reshape(p_GO_m.y,[],1)];
            Q_bs = -Q_bs;
        else
            b_eq = [b_eq; reshape(y_p_GO_m(:,1:size(y_m,2))-y_m,[],1)];
        end
        if n_c > 0
            Q_c = [Q_c; sparse(-repmat(Q_cs, n_s, 1))];
        end
        Q_b = blkdiag(Q_b, sparse(kron(eye(n_s),Q_bs))); % copy for each sample

        % inequality constraints Ax <= b
        R_a = [R_a; sparse(-repmat(R_as, n_s, 1))];
    end
end
R_c = sparse(size(R_a,1), n_c);
R_b = speye(size(R_a,1));

A_ineq = [R_a R_c R_b; R_a R_c -R_b];
b_ineq = sparse(size(A_ineq, 1),1);

A_eq = [sparse(size(Q_b,1), n_a) Q_c Q_b];

if strcmp(options_cs.task,"class") % classification
    % no equality constraints for classification constraints
    A_ineq = [A_eq; A_ineq];
    b_ineq = [b_eq; b_ineq];
    A_eq = []; %zeros(1,size(A_eq,2));
    b_eq = []; %0;
end

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
