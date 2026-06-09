function p_GO = computeGO(nlnsysDT,x0,u_ref,n_k, compute_params)
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
%    compute_params - boolean specifying if GO parameters are computed
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
% Written:       21-July-2023
% Last update:   12-February-2026 (LL, add input variable compute_params)
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

if nargin <= 4
    compute_params = true;
end

for k = 1 : n_k
    % compute reference solution x_ref and y_ref and the linearized system matrices
    x_ref(:,k+1) = nlnsysDT.mFile(x_ref(:,k), u_ref(:,k));
    y_ref(:,k) = nlnsysDT.out_mFile(x_ref(:,k), u_ref(:,k));

    % compute transfer matrices G for the x0->y(i) equation
    if compute_params
        % compute the linearized system matrices
        [A_lin{k},B_lin{k}] = nlnsysDT.jacobian(x_ref(:,k), u_ref(:,k));
        [C_lin_k,D_lin_k] = nlnsysDT.out_jacobian(x_ref(:,k), u_ref(:,k));

        p_GO = computeGO_k(A_lin, B_lin, C_lin_k, D_lin_k, k, p_GO);
    end
end
p_GO.x = x_ref;
p_GO.y = y_ref;
p_GO.u = u_ref;

% ------------------------------ END OF CODE ------------------------------
