function p_GO = computeGO_k(A_lin, B_lin, C_lin_k, D_lin_k, k, p_GO)
% computeGO_k - compute the GO parameters for a linearized state-space model
%       for time step k
%
% Syntax:
%    p_GO = computeGO_k(A_lin, B_lin, C_lin_k, D_lin_k, k, p_GO)
%
% Inputs:
%    A_lin - cell array with linearized A matrices for each time step
%    B_lin - cell array with linearized B matrices for each time step
%    C_lin - linearized C matrix for time step k
%    D_lin - linearized D matrix for time step k
%    k - current time step
%    p_GO - struct for saving the GO parameters
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
%    [1] L. Luetzow and M. Althoff, "Reachset-Conformant System
%        Identification," Transactions on Automatic Control, 2026. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       13-March-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

ny = size(C_lin_k,1);
nx = size(A_lin{1},1);

A_prod = eye(size(A_lin{1},1));
for j = 1 : k-1
    % iterate through past time steps
    A_prod = A_prod * A_lin{k-j};
    A_prod_j = 1;
    for i = 1 : k-j-1
        A_prod_j = A_prod_j * A_lin{k-i};
    end
    AA_prod = A_lin{k} * A_prod_j;
    p_GO.B{k,j} = AA_prod * B_lin{j};
    p_GO.F{k,j} = AA_prod * [eye(nx) zeros(nx, ny)];

    CA_prod = C_lin_k * A_prod_j;
    p_GO.D{k,j} = CA_prod * B_lin{j};
    p_GO.E{k,j} = CA_prod * [eye(nx) zeros(nx, ny)];
end

% save results for time step k
p_GO.A{k} = A_lin{k} * A_prod;
p_GO.B{k,k} = B_lin{k};
p_GO.F{k,k} = [eye(nx) zeros(nx, ny)]; % L = [L_x; L_y]

p_GO.C{k} = C_lin_k * A_prod;
p_GO.D{k,k} = D_lin_k;
p_GO.E{k,k} = [zeros(ny, nx) eye(ny)]; % L = [L_x; L_y]
end

% ------------------------------ END OF CODE ------------------------------
