function p_GO = computeGO(nlnARX,x0,u_ref,n_k)
% computeGO - compute the reference trajectory and the parameters for a 
%    linearized system
%
% Syntax:
%    p_GO = computeGO(nlnARX,x0,u_ref,n_k)
%
% Inputs:
%    nlnARX - system
%    x0 - stacked initial outputs
%    u_ref - reference input trajectory
%    n_k - number of time steps
%
% Outputs:
%    p_GO - struct with the GO parameters for a given reference trajectory
%               p_GO.A{k}      matrix that describes the influence of the 
%                              initial state x(1) on the state x(k+1)
%               p_GO.B{k,j}    matrix that describes the influence of the  
%                              input u(j) on the state x(k+1)     
%               p_GO.F{k,j}    matrix that describes the influence of the  
%                              linearization error L(j) on the state x(k+1)       
%               p_GO.C{k}      matrix that describes the influence of the 
%                              initial state x(1) on the output y(k)  
%               p_GO.D{k,j}    matrix that describes the influence of the  
%                              input u(j) on the output y(k)     
%               p_GO.E{k,j}    matrix that describes the influence of the  
%                              linearization error L(j) on the output y(k)
%               p_GO.x         reference state trajectory 
%                                   dimensions: n_x x (n_k+1)
%               p_GO.u         reference input trajectory  
%                                   dimensions: n_u x n_k
%               p_GO.y         reference output trajectory  
%                                   dimensions: n_y x n_k
%
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System
%        Identification," arXiv, 2024. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       14-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n_p = nlnARX.n_p;
n_y = nlnARX.nrOfOutputs;
n_u = nlnARX.nrOfInputs;
n_l = n_y;

if isa(x0, 'contSet')
    x0 = center(x0);
end
y_ref = zeros(n_y, n_k);
y_ref(:,1:n_p) = reshape(x0,n_y,[]);
u_stacked = getStackedU(nlnARX, u_ref);

% Linearization with y_d=y-y_ref and u_d=u-u_ref
%   x_d(k) = y_d(k-p+1:k)
%          = A_tilde{k} x_d(k-1) + \sum_{i=0}^{p} B_tilde{p+1-i,k} u_d(k-i)  
%               + F_tilde l(k)
%   y_d(k) = C_tilde x_d(k) 
A_tilde = cell(n_k,1);
B_tilde = cell(n_k,n_p+1);
F_tilde = [zeros(n_y*(n_p-1), n_l); eye(n_l)];
C_tilde = [zeros(n_y,(n_p-1)*n_y) eye(n_y)];
A_tilde_firstRows = [zeros(n_y*(n_p-1), n_y) eye(n_y * (n_p-1))];
B_tilde_firstRows = zeros(n_y*(n_p-1), n_u);

% Reformulation to
%   x_d(k) = A{k} x_d(p) + \sum_{i=0}^{k-1} B{i+1,k} u_d(i+1) 
%               + \sum_{i=0}^{k-1} F{i+1,k} l(k-i)
%   y_d(k) = C{k} x_d(p) + \sum_{i=0}^{k-1} D{i+1,k} u_d(i+1) 
%               + \sum_{i=0}^{k-1} E{i+1,k} l(k-i)
A = cell(n_k,1);
B = cell(n_k,n_k);
F = cell(n_k,n_k);
C = cell(n_k,1);
D = cell(n_k,n_k);
E = cell(n_k,n_k);

for k = n_p+1:n_k
    % compute linearized system matrices
    y_last = reshape(y_ref(:,k-n_p:k-1),[],1);
    [A_lin,B_lin] = nlnARX.jacobian(y_last, u_stacked(:,k));
    A_tilde{k} = [A_tilde_firstRows; A_lin];
    for i = 1:n_p+1
        B_tilde{k,i} = [B_tilde_firstRows; B_lin(:, (i-1)*n_u+1:i*n_u)];
    end
    y_ref(:,k) = nlnARX.mFile(y_last,  u_stacked(:,k));

    % compute reformulated system matrices
    A_k = cell(k-n_p,1);
    A_k{1} = A_tilde{k};
    F{k,k} = F_tilde;
    E{k,k} = C_tilde * F{k,k};
    for j = 1 : k-n_p-1
        A_k{j+1} = A_k{j} * A_tilde{k-j};
        F{k,k-j} = A_k{j} * F_tilde;
        E{k,k-j} = C_tilde * F{k,k-j};
    end
    A{k} = A_k{k-n_p};
    C{k} = C_tilde * A{k};

    % compute B and D
    for i=0:k-1
        if i <= n_p
            B{k,k-i} = B_tilde{k,n_p+1-i};
        else
            B{k,k-i} = zeros(size(B_tilde{k,n_p+1}));
        end

        for j = 1 : k-n_p-1
            if i >= j && i-j <= n_p
                B{k,k-i} = B{k,k-i} + A_k{j} * B_tilde{k-j,n_p+1-i+j};
            end
        end
        D{k,k-i} = C_tilde * B{k,k-i};
    end
end

% compute C and D
for k = 1:n_p
    C{k} = [zeros(n_y,(k-1)*n_y) eye(n_y) zeros(n_y,(n_p-k)*n_y)];
    [D{k,1:k}] = deal(zeros(n_y, n_u));
end

p_GO.x = zeros(n_p*n_y, n_k);
for i = 1:n_p
    p_GO.x((i-1)*n_y+1:i*n_y,n_p:end) = y_ref(:,i:end-n_p+i);
end

p_GO.y = y_ref;
p_GO.u = u_ref;

p_GO.A = A;
p_GO.B = B;
p_GO.F = F;

p_GO.C = C;
p_GO.D = D;
p_GO.E = E;

end

% ------------------------------ END OF CODE ------------------------------
