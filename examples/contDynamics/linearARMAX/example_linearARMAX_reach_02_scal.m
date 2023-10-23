function res = example_linearARMAX_reach_02_scal
% example_linearARMAX_reach_02_scal - example script to compare the
%   scalability of reachability algorithms for ARMAX models
%
% Syntax:
%    res = example_linearARMAX_reach_02_scal
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

% Authors:       Laura Luetzow
% Written:       04-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% user specifications
varying_param = "f_n"; % choose parameter which is varying
path_fig = ""; % path for saving the figure
num_rep = 5; % number of repetitions

% default parameter
fk_0 = 4;
fn_0 = 1;
p_0 = 2;

fk_vec = 1; 
fn_vec = 1; 
p_vec = 1;
if varying_param == "f_k"
    fk_vec = 1:1000:4001;
elseif varying_param == "f_n"
    fn_vec = 1:500:2001;
elseif varying_param == "p"
    p_vec = 1:75:301;
end

% time horizon
dt = 0.1;
params.tStart = 0;
rng('default')

% struct for saving the computation times
Ts = zeros(num_rep,1);
T = cell(1);
T{1,1} = 'exactAddition';
T{2,1} = 'tvpGeneral';
T{3,1} = 'tvpEfficient';
i_T = 2;

%% scalability analysis
for fk = fk_0*fk_vec
    for fn = fn_0*fn_vec
        dim_y = 2 * fn;
        dim_u = 3 * fn;

        A_bar = cell(1,1);
        B_bar = cell(1,1);
        B_bar{1} = randn(dim_y,dim_u);

        for p = p_0*p_vec
            kh = fk * p;

            for i = size(B_bar,1):p
                A_bar{i,1} = 0.01*randn(dim_y,dim_y);
                B_bar{i+1,1} = 0.01*randn(dim_y,dim_u);
            end

            % initialize linearSys objects
            sys_ARMAX = linearARMAX(A_bar, B_bar, dt);

            %% model parameters  ------------------------------------------

            % general parameters
            % disturbance set
            params.U = zonotope(zonotope(0.02+zeros(dim_u,1),0.01*diag(ones(dim_u,1))));            
            steps = kh + p;

            % input vector
            params.u = 0.01*randn(dim_u,steps+1);

            % noisy measurements which are used for the estimation
            params.y0 = rand(dim_y, p);
            params.tFinal = dt * steps;


            %% reach ------------------------------------------------------

            for i_m = 1 : 3    
                for i_j = 1:num_rep
                    tic;
                    options.armaxAlg = T{i_m,1};
                    R_ARMAX = reach(sys_ARMAX,params,options);
                    Ts(i_j) = toc;
                end
                T{i_m,i_T} = median(Ts);
            end
            T{4,i_T} = sprintf("k%d", kh);
            T{5,i_T} = sprintf("dim%d", fn);
            T{6,i_T} = sprintf("p%d", p);
            i_T = i_T + 1;
        end
    end
end

if varying_param == "f_k"
    vec = fk_0*fk_vec;
elseif varying_param == "f_n"
    vec = fn_0*fn_vec;
elseif varying_param == "p"
    vec = p_0*p_vec;
end
aux_plotTimes(vec, T, varying_param, path_fig)

% example completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_plotTimes(vec, T, varying_param, path_fig)
% function for plotting the computation times

col = colormap("lines");
lines = ["-", "--", ":", "-."];

figure('Position', [488 438 450 150]);
hold on; grid on;
box on
fontsize(gca,"scale", 1.1);
fontname(gca,"times");
axis = gca;
axis.LineWidth = 1;
plot(vec, cell2mat(T(1,2:end)), 'DisplayName', 'Proposition 3','LineWidth',1, ...
    'LineStyle', lines(1), 'Color', col(4,:), 'Marker', '*');
plot(vec, cell2mat(T(2,2:end)), 'DisplayName', 'Theorem 2','LineWidth',1, ...
    'LineStyle', lines(2), 'Color', col(5,:), 'Marker', '*');
plot(vec, cell2mat(T(3,2:end)), 'DisplayName', 'Theorem 3','LineWidth',1.5, ...
    'LineStyle', lines(3), 'Color', col(6,:), 'Marker', '*');
xlabel(varying_param);
ylabel("Time [s]");
lg = legend('Location', 'northwest');
if ~(isempty(path_fig) || path_fig == "")
    exportgraphics(gca,path_fig+sprintf('compTimes_%s.pdf',varying_param),'ContentType','vector');
end
end

% ------------------------------ END OF CODE ------------------------------
