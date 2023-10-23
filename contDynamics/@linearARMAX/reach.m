function [R,res] = reach(sys,params,options,varargin)
% reach - computes the reachable set for linear ARMAX models 
%
% Syntax:
%    R = reach(obj,reach_alg,options)
%
% Inputs:
%    sys - linearARMAX system
%    params - model parameters
%    options - options for the computation of reachable sets
%       options.armaxAlg - algorithm for reachability analysis 
%                   ('exactAddition', 'tvpEfficient' or 'tvpGeneral')
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       09-February-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% preprocessing of params
options = validateOptions(sys,mfilename,params,options);

% time period and number of steps
p = sys.dim;
tVec = options.tStart:sys.dt:options.tFinal;
n_y = sys.nrOfOutputs;

% initialize
Y = cell(length(tVec),1);
for i=1:p
    Y{i} = params.y0(:, i);   
end
U_const = params.U;
if size(params.u,2) == 1
    u = repmat(params.u, 1, length(tVec));
else 
    u = params.u;
end

% reachability analysis
if strcmp(options.armaxAlg, 'exactAddition') 
    % use exact addition for dependency-preserving sets [1, Proposiiton 3]

    U = cell(length(tVec),1);
    if isa(U_const, 'zonotope')
        % do linear transformation and exact addition directly on the
        % center vector and generator matrix and create zonotope at the end
        % (for higher efficiency)
        G_Uconst = generators(U_const);
        c_Uconst = u + center(U_const);
        G_Yfinal = cell(length(tVec),1);
        c_Yfinal = cell(length(tVec),1);
        for i=0:p-1
            G_Yfinal{i+1} = [];
            c_Yfinal{i+1} = Y{i+1};
        end
        eta = size(params.U.G,2);
        for k=p:length(tVec)-1
            %Y{k+1} = zonotope(zeros(n_y,(k+1)*eta+1));
            c_Y = 0;
            G_Y = zeros(n_y,(k+1)*eta);
            for i=1:p
                %Y_Y = zonotope(zeros(n_y,i*eta+1)) + sys.A_bar{i} * Y{k+1-i};
                %Y{k+1} = exactPlus(Y{k+1}, Y_Y);
                if isa(Y{k+1-i}, 'zonotope')
                    G_Y = G_Y + [zeros(n_y,i*eta) sys.A_bar{i} * G_Yfinal{k+1-i}];
                    c_Y = c_Y + sys.A_bar{i} * c_Yfinal{k+1-i};
                else
                    c_Y = c_Y + sys.A_bar{i} * Y{k+1-i};
                end
            end
            %Y_U = 0;
            c_U = 0;
            G_U = [];
            for i=0:p
                %Y_U = Y_U + sys.B_bar{i+1}*U{k-i+1};
                c_U = c_U + sys.B_bar{i+1}*c_Uconst(:,k-i+1);
                G_U = [G_U sys.B_bar{i+1}*G_Uconst];
            end
            %Y{k+1} = exactPlus(Y{k+1}, Y_U);
            c_Yfinal{k+1} = c_U + c_Y;
            G_Yfinal{k+1} = [G_U+G_Y(:,1:size(G_U,2)) G_Y(:, size(G_U,2)+1:end)];
            Y{k+1} = zonotope(c_Yfinal{k+1}, G_Yfinal{k+1});
        end

    elseif isa(U_const, 'polyZonotope')

        % create independent input sets
        id_start = 1;
        id_delta = length(U_const.id);
        for i=1:length(tVec)
            U{i} = replaceId(u(:, i) + U_const, id_start:id_start+id_delta-1);
            id_start = id_start + id_delta;
        end
        for i=1:p
            if ~isa(Y{i}, 'polyZonotope')
                Y{i} = polyZonotope(Y{i});
            end
        end

        % compute reachable set
        for k=p:length(tVec)-1
            Y{k+1} = sys.B_bar{1}*U{k+1};
            for i=1:p
                Y{k+1} = exactPlus(Y{k+1}, sys.A_bar{i} * Y{k+1-i});
                Y{k+1} = exactPlus(Y{k+1}, sys.B_bar{i+1} * U{k+1-i});
            end
        end
    else
        error('Other set representations not implemented yet');
    end

elseif strcmp(options.armaxAlg, 'tvpEfficient')
    % use efficient algorithm [1, Theorem 3]

    y_init = reshape(params.y0,[],1);
    u = [params.u zeros(size(params.u,1),p)];
    B_tilde = cell(2*p,1);
    if ~sys.tvp
        sys.setTVP;
    end

    % [1, Algorithm 2]:
    k = p;
    k_plus = k+p-1;
    computeTVP(sys,k,0);
    A_tilde = sys.A_tilde{k+1}; % [1, Eq. (9a)]
    B_tilde(1) = sys.B_tilde(1, k+1); % [1, Eq. (9b)]
    B_tilde = aux_compute_Btilde_Eq16(sys, B_tilde, 1:k_plus, k); % [1, Eq. (16)]

    % initialize
    S_c = 0;
    for i=0:k-1
        S_c = S_c + B_tilde{i+1} * U_const;
    end

    T_c1 = 0;
    for i=k:k_plus
        T_c1 = T_c1 + B_tilde{i+1} * U_const;
    end

    s_v = A_tilde * y_init;
    for i=p:k_plus
        s_v = s_v + B_tilde{i+1} * u(:,k_plus-i+1);
    end

    while true
        % compute Y_tilde
        y_v = s_v;
        for i=0:p-1
            y_v = y_v + B_tilde{i+1} * u(:,k_plus-i+1);
        end
        Y_tilde = S_c + T_c1 + y_v;

        % split Y_tilde in the output sets
        for j=0:p-1
            Y{k+1+j} = project(Y_tilde, j*n_y+1:(j+1)*n_y);
            if k+j >= length(tVec)-1
                break
            end
        end

        % check if final time point is reached
        if k_plus >= length(tVec)-1
            break
        end

        k = k + p;
        k_plus = k_plus + p;

        % compute new parameters 
        if k < 3*p
            computeTVP(sys,k,p);
            B_tilde(p+1) = sys.B_tilde(p+1, k+1);
            B_tilde = aux_compute_Btilde_Eq16(sys, B_tilde, p+1:2*p-1, k); % [1, Eq. (16)]
        end

        % initialize T_c2
        if k == 2*p
            T_c2 = 0;
            for i=k-p:k_plus-p
                T_c2 = T_c2 + B_tilde{i+1} * U_const;
            end
        end

        % update
        S_c = S_c + T_c2;
        T_c1 = A_tilde * T_c1;
        T_c2 = A_tilde * T_c2;
        s_v = A_tilde * s_v;
        for i=p:2*p-1
            s_v = s_v + B_tilde{i+1} * u(:,k_plus-i+1);
        end
    end

else 
    % use general algorithm [1, Theorem 2]

    y_init = reshape(params.y0,[],1);
    u = [params.u zeros(size(params.u,1),p)];
    B_tilde = cell(2*p,1);
    if ~sys.tvp
        sys.setTVP;
    end

    % [1, Algorithm 1]:
    k = p;
    k_plus = k+p-1;
    computeTVP(sys,k,0);
    A_tilde = sys.A_tilde{k+1}; % [1, Eq. (9a)]
    B_tilde(1) = sys.B_tilde(1, k+1); % [1, Eq. (9b)]
    B_tilde = aux_compute_Btilde_Eq16(sys, B_tilde, 1:k_plus, k); % [1, Eq. (16)]

    % initialize
    S = A_tilde * y_init;
    for i=p:k_plus
        S = S + B_tilde{i+1} * (u(:,k_plus-i+1) + U_const);
    end

    while true
        % compute Y_tilde
        Y_tilde = S;
        for i=0:p-1
            Y_tilde = Y_tilde + B_tilde{i+1} * (u(:,k_plus-i+1) + U_const);
        end

        % split Y_tilde in the output sets
        for j=0:p-1            
            Y{k+1+j} = project(Y_tilde, j*n_y+1:(j+1)*n_y);
            if k+j >= length(tVec)-1
                break
            end
        end

        % check if final time point is reached
        if k_plus >= length(tVec)-1
            break
        end

        k = k + p;
        k_plus = k_plus + p;

        % compute new parameters 
        if k < 3*p
            computeTVP(sys,k,p);
            B_tilde(p+1) = sys.B_tilde(p+1, k+1);
            B_tilde = aux_compute_Btilde_Eq16(sys, B_tilde, p+1:2*p-1, k); % [1, Eq. (16)]
        end

        % update S
        S = A_tilde * S;
        for i=p:2*p-1
            S = S + B_tilde{i+1} * (u(:,k_plus-i+1) + U_const);
        end
    end

end

% construct reachable set object
timePoint.set = Y;
timePoint.time = num2cell(tVec');
R = reachSet(timePoint);
end


% Auxiliary functions -----------------------------------------------------

function B_tilde_k = aux_compute_Btilde_Eq16(sys, B_tilde_k, i_vec, k)
p = sys.dim;
if isempty(sys.A_tilde{k+1})
    sys.A_tilde{k+1} = sys.conv_tvp.A_ext^(k); 
end

for i = i_vec
    B_tilde_k{i+1} = sys.conv_tvp.A_ext*B_tilde_k{i};
    if i+1 <= p+1
        B_tilde_k{i+1} = B_tilde_k{i+1} + sys.conv_tvp.B_ext{i+1};
    end
    if i+1-k > 0 && i+1-k <= p+1
        B_tilde_k{i+1} = B_tilde_k{i+1} - sys.A_tilde{k+1} * sys.conv_tvp.B_ext{i+1-k};
    end
end
end

% ------------------------------ END OF CODE ------------------------------
