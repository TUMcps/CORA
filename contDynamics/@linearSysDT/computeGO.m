function p_GO = computeGO(linsysDT,x0,u_ref,n_k)
% computeGO - compute the parameters of a general output (GO) model
%
% Syntax:
%    p_GO = computeGO(linsysDT,x0,u_ref,n_k)
%
% Inputs:
%    linsysDT - system
%    x0 - initial state
%    u_ref - reference input trajectory
%    n_k - number of time steps
%
% Outputs:
%    p_GO - struct with the GO parameters for a given reference trajectory
%               p_GO.A{k}      matrix that describes the influence of the 
%                               initial state x(1) on the state x(k+1)
%               p_GO.B{k,j}    matrix that describes the influence of the  
%                               input u(j) on the state x(k+1)     
%               p_GO.F{k,j}    matrix that describes the influence of the  
%                               linearization error L(j) on the state x(k+1)       
%               p_GO.C{k}      matrix that describes the influence of the 
%                               initial state x(1) on the output y(k)  
%               p_GO.D{k,j}    matrix that describes the influence of the  
%                               input u(j) on the output y(k)     
%               p_GO.E{k,j}    matrix that describes the influence of the  
%                               linearization error L(j) on the output y(k)
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

if ~isempty(x0) && ~isempty(u_ref)
    n_s = size(u_ref,3);
    x_ref = zeros(linsysDT.nrOfDims, n_k+1,n_s);
    y_ref = zeros(linsysDT.nrOfOutputs, n_k,n_s);
    if isa(x0, 'contSet')
        x0 = center(x0);
    end
    x_ref(:,1,:) = x0;
else
    x_ref = [];
    y_ref = [];
end
if ~isempty(u_ref)
    n_k = min(size(u_ref,2),n_k);
end
p_GO.A = cell(n_k,1);
p_GO.B = cell(n_k,n_k-1);
p_GO.F = cell(n_k,n_k-1);
p_GO.C = cell(n_k,1);
p_GO.D = cell(n_k,n_k);
p_GO.E = cell(n_k,n_k);

% compute the linearized system matrices
if iscell(linsysDT.A)
    %linear time-varying system
    for k = 1 : n_k

        if ~isempty(x0) && ~isempty(u_ref)
            % compute reference solution x_ref and y_ref
            x_ref(:,k+1,:) = pagemtimes(linsysDT.A{k},x_ref(:,k,:)) + pagemtimes(linsysDT.B{k},u_ref(:,k,:));
            y_ref(:,k,:) = pagemtimes(linsysDT.C{k},x_ref(:,k,:)) + pagemtimes(linsysDT.D{k},u_ref(:,k,:));
        end

        % compute transfer matrices G for the x0->y(i) equation
        A_prod = eye(size(linsysDT.A{1},1));
        for j = 1 : k-1
            A_prod = A_prod * linsysDT.A{k-j};
            A_prod_j = 1;
            for i = 1 : k-j-1
                A_prod_j = A_prod_j * linsysDT.A{k-i};
            end
            AA_prod = linsysDT.A{k} * A_prod_j;
            p_GO.B{k,j} = AA_prod * linsysDT.B{j};
            p_GO.F{k,j} = AA_prod * [eye(linsysDT.nrOfDims) zeros(linsysDT.nrOfDims, linsysDT.nrOfOutputs)];

            CA_prod = linsysDT.C{k} * A_prod_j;
            p_GO.D{k,j} = CA_prod * linsysDT.B{j};
            p_GO.E{k,j} = CA_prod * [eye(linsysDT.nrOfDims) zeros(linsysDT.nrOfDims, linsysDT.nrOfOutputs)];
        end
        p_GO.A{k} = linsysDT.A{k} * A_prod;
        p_GO.B{k,k} = linsysDT.B{k};
        p_GO.F{k,k} = [eye(linsysDT.nrOfDims) zeros(linsysDT.nrOfDims, linsysDT.nrOfOutputs)]; % L = [L_x; L_y]

        p_GO.C{k} = linsysDT.C{k} * A_prod;
        p_GO.D{k,k} = linsysDT.D{k};
        p_GO.E{k,k} = [zeros(linsysDT.nrOfOutputs, linsysDT.nrOfDims) eye(linsysDT.nrOfOutputs)]; % L = [L_x; L_y]
    end
else
    %linear time-invariant system    
    for k = 1 : n_k

        if ~isempty(x0) && ~isempty(u_ref)
            % compute reference solution x_ref and y_ref
            x_ref(:,k+1,:) = pagemtimes(linsysDT.A,x_ref(:,k,:)) + pagemtimes(linsysDT.B,u_ref(:,k,:));
            y_ref(:,k,:) = pagemtimes(linsysDT.C,x_ref(:,k,:)) + pagemtimes(linsysDT.D,u_ref(:,k,:));
        end

        % alternative computation
        if k>1
            p_GO.A{k} = linsysDT.A * p_GO.A{k-1};
            p_GO.C{k} = linsysDT.C * p_GO.A{k-1};
        else
            p_GO.A{k} = linsysDT.A;
            p_GO.C{k} = linsysDT.C;
        end
        p_GO.B{k,k} = linsysDT.B;
        p_GO.D{k,k} = linsysDT.D;
        for j = k-1 : -1: 1
            p_GO.B{k,j} = p_GO.A{k-j} * linsysDT.B;
            p_GO.D{k,j} = linsysDT.C* p_GO.B{k,j+1};
        end

    end
end

p_GO.x = x_ref;
p_GO.y = y_ref;
p_GO.u = u_ref;

% ------------------------------ END OF CODE ------------------------------
