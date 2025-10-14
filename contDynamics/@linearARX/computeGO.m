function p_GO = computeGO(linARX, x0, u_ref, n_k)
% computeGO - compute the reference trajectory and the parameters of a GO
%   model
%
% Syntax:
%    p_GO = computeGO(linARX, x0, u_ref, n_k)
%
% Inputs:
%    linARX - system
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
% Written:       27-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
n_p = linARX.n_p;
n_y = linARX.nrOfOutputs;
n_u = linARX.nrOfInputs;

if isa(x0, 'contSet')
    x0 = center(x0);
end

y_ref = zeros(n_y, n_k, size(x0,3));
y_ref(:,1:n_p,:) = reshape(x0,n_y,[], size(x0,3));

% A
A_tilde_firstRows = [zeros(n_y*(n_p-1), n_y) eye(n_y * (n_p-1))];
A_ext = zeros(n_y, n_p*n_y);
for i=0:n_p-1
    A_ext(:,i*n_y+1:(i+1)*n_y) = linARX.A_bar{n_p-i};
end
A_ext = [A_tilde_firstRows; A_ext];

% B
B_ext = cell(n_p+1);
B_tilde_firstRows = zeros(n_y*(n_p-1), n_u);
for i = 1:n_p+1
    B_ext{i} = [B_tilde_firstRows; linARX.B_bar{i}];
end

% init matrices
A = cell(n_k,1);
B = cell(n_k,n_k);
C = cell(n_k,1);
D = cell(n_k,n_k);
E = [zeros(n_y,(n_p-1)*n_y) eye(n_y)];

for k = n_p:n_k-1

    if ~isempty(x0) && ~isempty(u_ref)
        % compute reference solution x_ref and y_ref
        y_ref(:,k+1,:) = pagemtimes(linARX.B_bar{1}, u_ref(:, k+1,:));
        for i=1:n_p
            y_ref(:,k+1,:) = y_ref(:,k+1,:) + pagemtimes(linARX.A_bar{i,1}, y_ref(:, k-i+1,:)) ...
                + pagemtimes(linARX.B_bar{i+1}, u_ref(:, k-i+1,:));
        end
    end

    % compute reformulated system matrices
    A_ext_powerk = cell(k+1-n_p,1);
    A_ext_powerk{1} = 1;
    for j = 0 : k-n_p
        A_ext_powerk{j+2} = A_ext_powerk{j+1} * A_ext;
    end
    A{k+1} = A_ext_powerk{k+1-n_p+1};
    C{k+1} = E * A{k+1};

    for i=0:k
        % j = 0:
        if k-i <= n_p
            B{k+1,i+1} = B_ext{k-i+1};
        else
            B{k+1,i+1} = zeros(size(B_ext{n_p+1}));
        end
        
        % j > 0:
        for j = 1 : k-n_p
            if k-i-j >= 0 && k-i-j <= n_p
                B{k+1,i+1} = B{k+1,i+1} + A_ext_powerk{j+1} * B_ext{k-i-j+1};
            end
        end
        D{k+1,i+1} = E * B{k+1,i+1};
    end
end

% initial time steps
for k = 0:n_p-1
    C{k+1} = [zeros(n_y,k*n_y) eye(n_y) zeros(n_y,(n_p-k-1)*n_y)];
    [D{k+1,1:k+1}] = deal(zeros(n_y, n_u));
end

% save nominal signals in p_GO
p_GO.x = zeros(n_p*n_y, n_k);
for i = 1:n_p
    p_GO.x((i-1)*n_y+1:i*n_y,n_p:end) = y_ref(:,i:end-n_p+i);
end
p_GO.y = y_ref;
p_GO.u = u_ref;

% save matrices of GO model in p_GO
p_GO.A = A;
p_GO.B = B;
p_GO.F = 0;

p_GO.C = C;
p_GO.D = D;
p_GO.E = 0;

end

% ------------------------------ END OF CODE ------------------------------
