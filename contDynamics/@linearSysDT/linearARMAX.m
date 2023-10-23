function sys = linearARMAX(sys, method)
% linearARMAX - convert linear discrete-time system to linear ARMAX system
%
% Syntax:
%    sys = linearARMAX(sys,dt)
%
% Description:
%    Converts a discrete-time linear system to a equivalent discrete-time 
%    linear ARMAX system (the process and measurement disturbances w and v
%    of the linear system will be combined with the input u to u_ARMAX =
%    [u; w; v])
%
% Inputs:
%    sys - discrete-time linear system (class: linearSysDT)
%    dt - sampling time
%
% Outputs:
%    sys - discrete-time linear ARMAX system (class: linearARMAX)
%
% Example:
%    % discrete-time system
%    A = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
%    B = [0 0;0 0;1 0;0 1];
%    C = [1 0 0 0;0 1 0 0];
%    dt = 0.1;
%
%    sys = linearSysDT(A,B,[],C,dt);
% 
%    % convert to ARMAX system
%    sys_ARMAX = linearARMAX(sys);
%
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

% See also: linearSysDT, linearARMAX

% Authors:       Laura Luetzow
% Written:       03-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to discrete-time
dt = sys.dt;
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
E = eye(size(C,1)); %sys.E;
F = eye(size(A,1)); %sys.F;

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
    B = [B F zeros(dim_x,dim_y)];
    D = [D zeros(dim_y,dim_x) E];
    dim_u = size(B,2);

    % transform to ARMAX model
    tol = 1e4*eps;
    A_poly = poly(A);
    A_poly(abs(A_poly)<tol) = 0;
    p = size(A_poly,2) - 1;

    A_bar = cell(p,1);
    B_bar = cell(p+1,1);

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
        error("System is not observable.\n")
    end

    syms M_sym [dim_x dim_y]

    for k=1:10
        M_sol = solve((A+M_sym*C)^k == zeros(dim_x), M_sym);
        M_names = fieldnames(M_sol);
        if (dim_x+dim_y>2 && all(size(M_sol.(M_names{1,1})) == [1,1])) || ...
                (dim_x+dim_y==2 && all(size(M_sol) == [1,1]))
            break
        end
    end
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
    if sum(abs((A+M*C)^k), 'all') > 1e-5
        warning("Low Accuracy of Transformation matrix M.")
    end

    % compute ARMAX parameters
    A_bar = cell(p,1);
    B_bar = cell(p+1,1);
    B_bar{1} = [D zeros(dim_y,dim_x) E];

    for k=1:p
        A_bar{k,1} = -C*(A+M*C)^(k-1)*M;
        B_bar{k+1,1} = C*(A+M*C)^(k-1)* [(B+M*D) F M*E];
    end
end
sys = linearARMAX(A_bar, B_bar, dt);
    
end

% ------------------------------ END OF CODE ------------------------------
