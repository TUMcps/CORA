function linARX = linearARX(linsysDT, method, E_x, E_y)
% linearARX - convert linear discrete-time system to linear ARX system
%
% Syntax:
%    linARX = linearARX(linsysDT,method)
%    linARX = linearARX(linsysDT,method,E)
%    linARX = linearARX(linsysDT,method,E,F)
%
% Description:
%    Converts a discrete-time linear system to a equivalent discrete-time 
%    linear ARX system (the process and measurement disturbances w and v
%    of the linear system will be combined with the input u to u_ARX =
%    [u; w; v])
%
% Inputs:
%    linsysDT - discrete-time linear system (class: linearSysDT)
%    method - conversion method
%    E_x - process noise matrix
%        default: identity
%    E_y - measurement noise matrix
%        default: identity
%
% Outputs:
%    linARX - discrete-time linear ARX system (class: linearARX)
%
% Example:
%    % discrete-time system
%    A = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
%    B = [0 0;0 0;1 0;0 1];
%    C = [1 0 0 0;0 1 0 0];
%    dt = 0.1;
%    linsysDT = linearSysDT(A,B,[],C,dt);
% 
%    % convert to ARX system
%    linARX = linearARX(linsysDT);
% 
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.
%   [2] M. Aoki and A. Havenner, "State space modeling of multiple time 
%       series," Econometric Reviews, vol. 10, no. 1, pp. 1-59, 1991.
%   [3] M. Phan, L. Horta, J.-N. Juang, and R. Longman, “Linear system
%       identification via an asymptotically stable observer,” Journal of
%       Optimization Theory and Application, vol. 79, no. 1, pp. 59–86, 
%       1993.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT, linearARX

% Authors:       Laura Luetzow
% Written:       03-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to discrete-time
dt = linsysDT.dt;
A = linsysDT.A;
B = linsysDT.B;
C = linsysDT.C;
D = linsysDT.D;
if nargin < 4
    E_y = eye(size(C,1)); %sys.E;
end
if nargin < 3
    E_x = eye(size(A,1)); %sys.F;
end

dim_x = size(C,2);
dim_y = size(C,1);

if isempty(D)
    D = zeros(dim_y, size(B,2));
end

% choose transformation method
if nargin ~= 2
    if dim_x > 4
        method = "CH";
    else
        method = "M";
    end
end

if method == "CH"
    % transformation method based on Cayley-Hamiliton theorem, see [2]
    % (does not find the minimal model order but scales well to high
    % dimensions)

    % extend input with disturbances u_new = [u; w; v]
    B = [B E_x zeros(dim_x,size(E_y,2))];
    D = [D zeros(dim_y,size(E_x,2)) E_y];
    dim_u = size(B,2);

    tol = 1e4*eps;

    % transform to ARX model
    A_poly = poly(A);
    A_poly(abs(A_poly)<tol) = 0;
    p = size(A_poly,2) - 1;

    % initialize ARX matrices A_bar and B_bar
    A_bar = cell(p,1);
    B_bar = cell(p+1,1);

    % compute A_bar and B_bar
    CABs = {D};
    for i = 1:p
        A_bar{i} = -A_poly(i+1)*eye(dim_y);
        CABs = [{C*A^(i-1)*B}; CABs];
    end

    for i = 0:p
        B_bar{i+1} = CABs{end-i};
        for j=1:i
            B_bar{i+1} = B_bar{i+1} + A_poly(j+1)*CABs{end-i+j};
        end
    end

else
    % transformation using the deadbeat observer gain M, see [3]
    % (finds minimal model order p but does not scale well for high
    % dimensions)

    %find transformation matrix M
    Ob = C;
    for i=1:dim_x-1
        Ob = [Ob; C*A^i];
    end
    if rank(Ob) ~= dim_x
        throw(CORAerror('CORA:notSupported',...
            "System is not observable.\n"));
    end

    % create symbolic matrix M_sym
    syms M_sym [dim_x dim_y]

    for k=1:10
        % solve for M_sym
        M_sol = solve((A+M_sym*C)^k == zeros(dim_x), M_sym);
        M_names = fieldnames(M_sol);
        if (dim_x+dim_y>2 && all(size(M_sol.(M_names{1,1})) == [1,1])) || ...
                (dim_x+dim_y==2 && all(size(M_sol) == [1,1]))
            break
        end
    end

    % compute M from M_sol
    M = zeros(dim_x, dim_y);
    if dim_x > 1
        for i_x=1:dim_x
            if dim_y > 1
                for i_y = 1:dim_y
                    M(i_x, i_y) = double(M_sol.(sprintf("M_sym%d_%d", i_x, i_y)));
                end
            else
                M(i_x, 1) = double(M_sol.(M_names{i_x, 1}));
            end
        end
    else
        M(1, 1) = double(M_sol);
    end
    p = k; 
    
    % check accuracy
    if sum(abs((A+M*C)^k), 'all') > 1e-5
        CORAwarning('CORA:contDynamics',"Low Accuracy of Transformation matrix M.")
    end

    % compute ARX parameters
    A_bar = cell(p,1);
    B_bar = cell(p+1,1);
    B_bar{1} = [D zeros(dim_y,size(E_x,2)) E_y];

    if ~isempty(E_y)
        ME = M*E_y;
    else
        ME = [];
    end
    for k=1:p
        A_bar{k,1} = -C*(A+M*C)^(k-1)*M;
        B_bar{k+1,1} = C*(A+M*C)^(k-1)* [(B+M*D) E_x ME];
    end
end

% instantiate linearARX object
linARX = linearARX(A_bar, B_bar, dt);

% ------------------------------ END OF CODE ------------------------------
