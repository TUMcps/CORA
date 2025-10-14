function [sys, R0, U, p_true] = loadDynamics(dynamics, type, p)
% loadDynamics - load system dynamics and uncertainty sets
%
% Syntax:
%    [sys, R0, U] = loadDynamics(dynamics)
%
% Inputs:
%    dynamics - string specifying the dynamical system
%    type - string "standard" (default) for the normal uncertainty sets or
%               "diag" for diagonal generator matrices with random elements
%               "rand" for non-diagonal generator matrices
%    p - [optional] model parameters 
%
% Outputs:
%    sys - dynamical system
%    R0 - default initial state set
%    U - default input set
%    p_true - true model parameters
%
% References:
%   [1] L. Luetzow, M. Althoff, "Reachability analysis of ARMAX models," in 
%       Proc. of the 62nd IEEE Conference on Decision and Control, pp. 
%       7021–7028, 2023.
%   [2] E. N. Lorenz, “Deterministic nonperiodic flow,” Journal of
%       Atmospheric Sciences, vol. 20, no. 2, pp. 130 – 141, 1963.
%   [3] A. Kroll and H. Schulte, “Benchmark problems for nonlinear system
%       identification and control using soft computing methods: Need and
%       overview," Applied Soft Computing, vol. 25, pp. 496–513, 2014.
%   [4] J. M. Bravo. "Robust MPC of constrained discrete-time nonlinear 
%       systems based on approximated reachable sets", Automatica, 2006.
%   [5] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization", in Proc. 
%       of the 62nd IEEE Conference on Decision and Control, pp.
%       4042-4048, 2008.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: testCase

% Authors:       Laura Luetzow
% Written:       01-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    type = "standard";
end

switch dynamics
    case "pedestrian"
        % pedestrian model as a state-space model [1]
        p_true = [1 0.01 5e-5 0.01]';
        if nargin < 3
            p = p_true;
        end
        A = [p(1)	0	    p(2)	0
            0	    p(1)	0	    p(2)
            0	    0	    p(1)	0
            0	    0	    0	    p(1)];
        B =[p(3)    0       0       0
            0	    p(3)    0       0
            p(4)    0       0       0
            0	    p(4)    0       0];
        C =[1	    0	    0	    0
            0	    1	    0	    0];
        D =[0	    0       1       0
            0	    0       0       1];

        dt = 0.01;
        sys = linearSysDT(A,B,[],C,D, dt);

        % construct uncertainty sets
        dim_x = length(sys.A);
        dim_u = 2;
        dim_v = 2;
        switch type
            case "rand" 
                % random uncertainty sets
                c_R0 = [-0.76; -9.68; 0.21; -5.42];
                c_U = [-0.16; -8.93];
                c_V = [1.48; -7.06];
                G_R0 = [-0.02 0.13 0.10 0.06
                    0.30 -0.24  0.21 -0.16
                    0.28  0.14  0.15  0.18
                    0.28  0.33 -0.06 -0.23];
                G_U = [0.07   -0.25
                    -0.28   -0.11];
                G_V = [-0.08    0.01
                    -0.00   -0.03];
            case "diag" 
                % diagonal uncertainty sets
                c_R0 = 0.1*[-0.76; -9.68; 0.21; -5.42];
                c_U = 0.1*[-0.16; -8.93];
                c_V = 0.1*[1.48; -7.06];
                G_R0 = diag([0.22 0.13 0.10 0.06]);
                G_U = diag([0.07 0.25]);
                G_V = diag([0.08 0.01]);
            case "standard" 
                % standard uncertainty sets
                c_R0 = zeros(dim_x,1);
                G_R0 = [];
                c_U = 0.1+zeros(dim_u,1);
                G_U = 0.2*diag(ones(dim_u,1));
                c_V =  -0.05+zeros(dim_v,1);
                G_V = 0.1*[diag(ones(dim_v,1)) ones(dim_v,1)];
        end
        R0 = zonotope([c_R0,G_R0]);
        V = zonotope([c_V,G_V]);
        U = cartProd(zonotope([c_U,G_U]), V);

    case "pedestrianARX"
        % pedestrian model as an ARX model [1]
        p_true = [2 -1 5e-5 -2]';
        if nargin < 3
            p = p_true;
        end
        A{1,1} = [  p(1)	0	    
                    0	    p(1)];
        A{2,1} = [  p(2)	0	    
                    0	    p(2)];
        B{1,1} = [  0	    0       1       0
                    0	    0       0       1];
        B{2,1} = [  p(3)    0       p(4)    0
                    0	    p(3)    0       p(4)];
        B{3,1} = [  p(3)    0       1       0
                    0	    p(3)    0       1];
        dt = 0.01;
        sys = linearARX(A, B, dt);

        % construct uncertainty sets
        dim_x = 4;
        dim_u = 2;
        dim_v = 2;
        switch type
            case "rand" 
                % random uncertainty sets
                c_U = [-0.16; -8.93];
                c_V = [1.48; -7.06];
                G_U = [0.07   -0.25
                    -0.28   -0.11];
                G_V = [-0.08    0.01
                    -0.00   -0.03];
            case "diag" 
                % diagonal uncertainty sets
                c_U = 0.1*[-0.16; -8.93];
                c_V = 0.1*[1.48; -7.06];
                G_U = diag([0.07 0.25]);
                G_V = diag([0.08 0.01]);
            case "standard" 
                % standard uncertainty sets
                c_U = 0.1+zeros(dim_u,1);
                G_U = 0.2*diag(ones(dim_u,1));
                c_V =  -0.05+zeros(dim_v,1);
                G_V = 0.1*[diag(ones(dim_v,1)) ones(dim_v,1)];
        end
        R0 = zonotope(zeros(dim_x,1));
        V = zonotope([c_V,G_V]);
        U = cartProd(zonotope([c_U,G_U]), V);            

    case "lorenz"
        % Lorenz system [2]
        p_true = [10 28 8/3]';
        if nargin < 3
            p = p_true;
        end
        dt = 0.01;
        fun = @(x,u) aux_dynLorenz(x,u,dt,p);
        dim_x = 3;
        dim_u = 3;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y);
        sys = nonlinearSysDT('lorenz', fun, dt, dim_x, dim_u, out_fun, dim_y);

        % construct uncertainty sets
        switch type 
            case "rand" 
                % random uncertainty sets
                c_R0 = [6.01; 9.36; -3.73];
                c_U = [7.56; -8.03; -1.57];
                G_R0 = [-0.11   -0.11    0.04
                    -0.07    0.08    0.14
                    0.02     0       0.02];
                G_U = [0.06   -0.02    0.04
                    -0.04    0.08   -0.09
                    -0.07    0.20   -0.10];
            case "diag" 
                % diagonal uncertainty sets
                c_R0 = 0.1*[6.01; 9.36; -3.73];
                c_U = 0.1*[7.56; -8.03; -1.57];
                G_R0 = 0.03*diag([0.11   0.11    0.24]);
                G_U = diag([0.06   -0.02    0.04]);
            case "standard" 
                % standard uncertainty sets
                c_U = [0.5;0.1;-0.2];
                c_R0 = [2; -1; 4];
                G_U =  diag([0.1;2;0.2]);
                G_R0 = 0.2*eye(dim_x);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_U,G_U]);
        V = zonotope([]);
        U = cartProd(W, V);
    
    case "lorenz_2D"
        % first two dimensions of the Lorenz system [2]
        p_true = [10 28]';
        if nargin < 3
            p = p_true;
        end
        dt = 0.01;
        fun = @(x,u) aux_dynLorenz2D(x,u,dt,p);
        dim_x = 2;
        dim_u = 2;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y);
        sys = nonlinearSysDT('lorenz', fun, dt, dim_x, dim_u, out_fun, dim_y);

        % construct uncertainty sets
        switch type 
            case "rand" 
                % random uncertainty sets
                c_R0 = [6.01; 9.36];
                c_U = [7.56; -8.03];
                G_R0 = [-0.11   -0.11
                    -0.07    0.08];
                G_U = [0.06   -0.02
                    -0.04    0.08];
            case "diag"
                % diagonal uncertainty sets
                c_R0 = 0.1*[6.01; 9.36];
                c_U = 0.1*[7.56; -8.03];
                G_R0 = 0.03*diag([0.11   0.11]);
                G_U = diag([0.06   -0.02]);
            case "standard" 
                % standard uncertainty sets
                c_U = [0.5;0.1];
                c_R0 = [2; -1];
                G_U =  diag([0.1;2]);
                G_R0 = 0.2*eye(dim_x);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_U,G_U]);
        V = zonotope([]);
        U = cartProd(W, V);

    case "NARX" 
        % artificial NARX model, adapted from [3]
        p_true = [0.8 1.2]';
        if nargin < 3
            p = p_true;
        end

        f = @(y,u) [y(1,1)/(1+y(2,1)^2) + p(1)*u(3,1); ...
            (y(1,1) * y(2,1))/(1+y(2,1)^2)+ p(2)*u(6,1)];
        dt = 0.1;
        dim_y = 2;
        dim_u = 2;
        p_dim = 2;
        sys = nonlinearARX(dynamics,f,dt,dim_y, dim_u, p_dim);

        % construct uncertainty sets
        R0 = zonotope(zeros(dim_y*p_dim,1));
        switch type
            case "rand"
                % random uncertainty sets
                c_U = [-1.66; 4.41];
                G_U = [-0.1   0.13
                    0.25   -0.09];
            case "diag" 
                % diagonal uncertainty sets
                c_U = 0.1*[-1.66; 4.41];
                G_U = 0.7*diag([0.1   0.13]);
            case "standard" 
                % standard uncertainty sets
                c_U = [0;0.05];
                G_U =  0.2*eye(2);
        end
        U = zonotope(c_U, G_U);

    case "Square" 
        % artificial, simple NARX model
        p_true = [0.8 1.2]';
        if nargin < 3
            p = p_true;
        end

        f = @(y,u) [y(1,1)^2 + p(1)*u(3,1); ...
            y(2,1)^2+p(2)*u(2,1)];
        dt = 0.1;
        dim_y = 2;
        dim_u = 2;
        p_dim = 1;
        sys = nonlinearARX(dynamics,f,dt,dim_y, dim_u, p_dim);

        % construct uncertainty sets
        R0 = zonotope(zeros(dim_y*p_dim,1));
        switch type
            case "rand" 
                % random uncertainty sets
                c_U = [-1.66; 4.41];
                G_U = [-0.1   0.13
                    0.25   -0.09];
            case "diag" 
                % diagonal uncertainty sets
                c_U = 0.1*[-1.66; 4.41];
                G_U = 0.7*diag([0.1   0.13]);
            case "standard" 
                % standard uncertainty sets
                c_U = [0;0.05];
                G_U =  0.2*eye(2);
        end
        U = zonotope(c_U, G_U);

    case "bicycle"
        % bicycle dynamics (see DOTBicycleDynamics_SRX_velEq.m)
        dt = 0.001;
        p_true = []; % no parameters defined
        fun = @(x,u) x + dt*DOTBicycleDynamics_SRX_velEq(x,u);
        dim_x = 6;
        dim_u = 4;
        out_fun = @(x,u) [x(4); x(5)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('bicycle', fun, dt, dim_x, dim_u, out_fun, dim_y);

        % construct uncertainty sets
        if nargin > 1 && type == "rand" 
            % random uncertainty sets
            c_R0 = [1.86; 3.46; 3.97; 5.39; 4.19; 6.85];
            c_W = [2.04; 8.78];
            c_V = [-0.02; 0.06];
            G_R0 = 0.01*[1.78 -1.37  0.79  0.60 -1.17 -1.4800
                1.77 -0.29  0.93 -0.54 -0.69  0.26
                -1.87  1.27 -0.49 -0.16  0.93 -2.02
                -1.05  0.07  1.80  0.61 -1.48  0.20
                -0.42  0.45  0.59 -1.04 -0.56  0.43
                1.40 -0.32 -0.64 -0.35 -0.03 -1.27];
            G_W = 0.01*diag([0.55 0.17]);
            G_W(1,2) = 0.2; 
            G_V = 0.002*eye(dim_y);
        else 
            % diagonal or standard uncertainty sets
            c_R0 = [1.2;0.5; 0; 0; 0; 0];
            c_W = zeros(dim_u-dim_y,1);
            c_V = zeros(dim_y,1);
            G_R0 = eye(6);
            G_W = eye(dim_u-dim_y);
            G_V = eye(dim_y);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "bicycleHO"
        % higher-order bicycle dynamics (see highorderBicycleDynamics.m)
        dt = 0.001;
        p_true = []; % no parameters defined
        fun = @(x,u) x + dt*highorderBicycleDynamics(x,u);
        dim_x = 18;
        dim_u = 4;
        out_fun = @(x,u) [x(5); x(6)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('bicycleHO', fun, dt, dim_x, dim_u, out_fun, dim_y);

        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([[1.2;0.5; 0; 5; zeros(14,1)],0.01*eye(dim_x)]);
        W = zonotope([zeros(2,1),0.004*eye(2)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

    case "cstrDiscr"
        % stirred-tank reactor system [5]        
        dt = 0.015;
        p_true = []; % no parameters defined
        fun = @(x,u) cstrDiscr(x,u,dt);
        dim_x = 2;
        dim_u = 4;
        out_fun = @(x,u) [x(1); x(2)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('cstrDiscr', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        % construct uncertainty sets
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        c_R0 = [-0.15;-45];
        c_W = zeros(dim_x,1);
        c_V = zeros(dim_y,1);
        G_R0 = diag([0.005;3]);
        G_W = diag([0.1;2]);
        G_V = 0.002*eye(dim_y);
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "tank"
        % stirred-tank reactor system with 6 dimensions [5]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank6EqDT(x,u,dt);
        dim_x = 6;
        dim_u = 4;
        out_fun = @(x,u) [x(1); x(2)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('tank6', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        % construct uncertainty sets
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        c_R0 = [2; 4; 4; 2; 10; 4];
        c_W = zeros(2,1);
        c_V = zeros(dim_y,1);
        G_R0 = 0.2*eye(6);
        G_W = diag([0.1;2]);
        G_V = 0.002*eye(dim_y);
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "tank30"
        % stirred-tank reactor system with 30 dimensions [5]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank30EqDT_inflow15(x,u,dt);
        dim_x = 30;
        dim_u = 15;
        dim_y = 6;
        out_fun = @(x,u) x(1:dim_y) + u(1:dim_y);
        sys = nonlinearSysDT('tank30', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        % construct uncertainty sets
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([12*rand(dim_x,1),0.2*eye(dim_x)]);
        W = zonotope([zeros(dim_y-dim_y,1),0.01*eye(dim_u-dim_y)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

    case "tank60"
        % stirred-tank reactor system with 60 dimensions [4]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank60EqDT_inflow30(x,u,dt);
        dim_x = 60;
        dim_u = 30;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y) + u(1:dim_y);
        sys = nonlinearSysDT('tank60', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        % construct uncertainty sets
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([12*rand(dim_x,1),0.2*eye(dim_x)]);
        W = zonotope([zeros(dim_y-dim_y,1),0.01*eye(dim_u-dim_y)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

end
end


% Auxiliary functions -----------------------------------------------------

function xnew = aux_dynLorenz2D(x,u,dt,p)
% dynamics of the Lorenz system

xdot = [(p(1)+u(1))*(x(2)-x(1)); ...
    (p(2)+u(2))*x(1)-x(2)-x(1)];
xnew = x + dt*xdot; 
end

function xnew = aux_dynLorenz(x,u,dt,p)
% dynamics of the Lorenz system

xdot = [(p(1)+u(1))*(x(2)-x(1)); ...
    (p(2)+u(2))*x(1)-x(2)-x(1)*x(3); ...
    x(1)*x(2)-(p(3)+u(3))*x(3)];
xnew = x + dt*xdot; 
end

% ------------------------------ END OF CODE ------------------------------
