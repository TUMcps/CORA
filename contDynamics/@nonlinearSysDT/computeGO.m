function p_GO = computeGO(nlnsysDT,x0,u_ref,n_k)
% computeGO - compute the reference trajectory and the parameters for a 
%    linearized system
%
% Syntax:
%    p_GO = computeGO(nlnsysDT,x0,u_ref,n_k)
%
% Inputs:
%    nlnsysDT - system
%    x0 - initial state
%    u_ref - reference input trajectory
%    n_k - number of time steps
%
% Outputs:
%    p_GO - struct with the GO parameters for a give nreference trajectory
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
% Written:       21-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x_ref = zeros(nlnsysDT.nrOfDims, n_k+1);
y_ref = zeros(nlnsysDT.nrOfOutputs, n_k);
if isa(x0, 'contSet')
    x0 = center(x0);
end
x_ref(:,1) = x0;

A_lin = cell(n_k-1,1);
B_lin = cell(n_k-1,1);

p_GO.A = cell(n_k,1);
p_GO.B = cell(n_k,n_k-1);
p_GO.F = cell(n_k,n_k-1);
p_GO.C = cell(n_k,1);
p_GO.D = cell(n_k,n_k);
p_GO.E = cell(n_k,n_k);

% compute the linearized system matrices
for k = 1 : n_k
    % compute reference solution x_ref and y_ref and the linearized system matrices
    [A_lin{k},B_lin{k}] = nlnsysDT.jacobian(x_ref(:,k), u_ref(:,k));
    x_ref(:,k+1) = nlnsysDT.mFile(x_ref(:,k), u_ref(:,k));

    [C_lin,D_lin] = nlnsysDT.out_jacobian(x_ref(:,k), u_ref(:,k));
    y_ref(:,k) = nlnsysDT.out_mFile(x_ref(:,k), u_ref(:,k));

    % compute transfer matrices G for the x0->y(i) equation
    A_prod = eye(size(A_lin{1},1));
    for j = 1 : k-1
        A_prod = A_prod * A_lin{k-j};
        A_prod_j = 1;
        for i = 1 : k-j-1
            A_prod_j = A_prod_j * A_lin{k-i};
        end
        AA_prod = A_lin{k} * A_prod_j;
        p_GO.B{k,j} = AA_prod * B_lin{j};
        p_GO.F{k,j} = AA_prod * [eye(nlnsysDT.nrOfDims) zeros(nlnsysDT.nrOfDims, nlnsysDT.nrOfOutputs)];

        CA_prod = C_lin * A_prod_j;
        p_GO.D{k,j} = CA_prod * B_lin{j};
        p_GO.E{k,j} = CA_prod * [eye(nlnsysDT.nrOfDims) zeros(nlnsysDT.nrOfDims, nlnsysDT.nrOfOutputs)];
    end
    p_GO.A{k} = A_lin{k} * A_prod;
    p_GO.B{k,k} = B_lin{k}; 
    p_GO.F{k,k} = [eye(nlnsysDT.nrOfDims) zeros(nlnsysDT.nrOfDims, nlnsysDT.nrOfOutputs)]; % L = [L_x; L_y]
    
    p_GO.C{k} = C_lin * A_prod;
    p_GO.D{k,k} = D_lin;
    p_GO.E{k,k} = [zeros(nlnsysDT.nrOfOutputs, nlnsysDT.nrOfDims) eye(nlnsysDT.nrOfOutputs)]; % L = [L_x; L_y]

end
p_GO.x = x_ref;
p_GO.y = y_ref;
p_GO.u = u_ref;

% ------------------------------ END OF CODE ------------------------------
