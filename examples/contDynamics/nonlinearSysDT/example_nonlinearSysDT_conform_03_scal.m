function completed = example_nonlinearSysDT_conform_03_scal
% example_nonlinearSysDT_conform_03_scal - - example for conformance 
%   identification of nonlinear discrete-time systems to analyze the 
%   scalability
% 
%
% Syntax:
%    example_nonlinearSysDT_conform_03_scal
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] 

% Authors:       Laura Luetzow
% Written:       21-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% User Specifications ----------------------------------------------------
%
% modus - Identification modus to choose:
%   - "N_m": vary the number of nominal trajectories
%   - "N_s": vary the number of samples per nominal trajectory
%   - "N_k": vary the length of the time horizon
%   - "N_n": vary the system dimension
% num_rep - number of repetitions for each specification

modus = "N_k";
num_rep = 10;

% Set Random Number Stream
rng('default')

% Values for the Constant Parameters
N_m_vec = 10; 
N_s_vec = 1; 
N_k_vec = 4; 
N_n_vec = 3;

% Values for the Varying Parameters
switch modus
    case "N_m"
        N_m_vec = [10 100 1000];
    case "N_s"
        N_s_vec = [1 10 100];
    case "N_k"
        N_k_vec = [4 20 40];
    case "N_n"
        N_n_vec = [3 6 9 100];
end

% Reachability Settings  
options_reach.zonotopeOrder = 100;
options_reach.tensorOrder = 2;
options_reach.errorOrder = 5;
options_reach.tensorOrderOutput = 2;

% Conformance Settings
options = options_reach;
options.cs.cost = 'interval';
options.cs.robustnessMargin = 1e-9;
constraints = {'half','gen'};

% Setup Table for Saving the Results
varTypes = {'double','double','double','double','double','double',...
    'double','double','double','double'};
varNames = ["N_m", "N_s", "N_k", "N_n", constraints{1}+"_mean", ...
    constraints{2}+"_mean", constraints{1}+"_std", ...
    constraints{2}+"_std", constraints{1}+"_med", constraints{2}+"_med"];
tTab = table('Size',[3 10],'VariableTypes', varTypes,'VariableNames',...
    varNames);
Ts = inf*ones(length(constraints),num_rep);
i_tab = 1;
maxTimeStand = 600;

for N_m = N_m_vec
    for N_s = N_s_vec
        for N_k = N_k_vec
            for N_n = N_n_vec

                % Parameters and System Dynamics 
                [sys, params_true.U, params_true.R0] = ...
                    aux_load_tank(N_n);
                params_true.tFinal = sys.dt * N_k - sys.dt;

                % Simulation 
                params_true.testSuite = aux_create_testSuite(sys, ...
                    params_true, N_m, N_s, N_k);

                % Initial Estimates of the Disturbance Sets
                params_id_init = params_true;
                params_id_init.R0 = zonotope([center(params_true.R0) ...
                    eye(sys.nrOfStates) ones(sys.nrOfStates,1)]);
                params_id_init.U = zonotope([center(params_true.U) ...
                    eye(sys.nrOfInputs) ones(sys.nrOfInputs,1)]);

                %% Conformance Identification -----------------------------

                for i_c = 1:length(constraints)
                    options.cs.constraints = constraints{i_c};

                    fprintf("Identification with %s-constraints, " + ...
                        "Nm=%d, Ns=%d, Nk=%d, Nx=%d\n",...
                        options.cs.constraints, N_m, N_s, N_k, N_n);
                    for j = 1:num_rep
                        f = parfeval(@conform,2,sys,params_id_init,options);
                        while f.State == "queued"
                        end
                        tic;
                        while f.State == "running"
                            if toc > maxTimeStand
                                cancel(f)
                            end
                        end
                        if f.State == "finished"
                            Ts(i_c,j) = toc;
                        end
                        if Ts(i_c,j) > maxTimeStand
                            Ts(i_c,:) = Inf;
                            break
                        end
                    end
                    fprintf("Median identification time: %.4f\n", ...
                        median(Ts(i_c,:)));
                end

                % Save the Computation Time
                tTab(i_tab,:) = {N_m, N_s, N_k, N_n, mean(Ts(1,:)), ...
                    mean(Ts(2,:)), std(Ts(1,:)), std(Ts(2,:)), ...
                    median(Ts(1,:)), median(Ts(2,:))};
                i_tab = i_tab + 1;
            end
        end
    end
end

% example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys, U, R0] = aux_load_tank(N_n)
% load the tank dynamics with the specified dimension

dynamics = sprintf('tank%d', N_n);
dt = 0.5;
N_u = ceil(N_n/3);
N_y = N_n;
fun = @(x,u) tankN (x,u,dt,N_n);
out_fun = @(x,u) x(1:N_y);
sys = nonlinearSysDT(dynamics, fun, dt, N_n, N_u, out_fun, N_y);
R0 = zonotope([10*rand(N_n,1),0.2*eye(N_n)]);
U = zonotope([100*rand(N_u,1),0.1*eye(N_u)]);
end

function testSuite = aux_create_testSuite(sys, params, N_m, N_s, N_k)
% create identification data by simulation (nominal values u_m and x0_m
% must be positive for tank system)

testSuite = cell(N_m,1);    
for m = 1:N_m
    y_m = cell(N_s,1);
    % sample a positive random nominal input trajectory 
    u_m = 10*rand(sys.nrOfInputs, N_k); 
    % sample a positive random nominal initial state vector
    x0_m = 10*rand(size(params.R0.center));

    % sample the true input trajectory and initial state vector
    for s = 1:N_s
        if rand(1) < 0.3
            params.x0 = x0_m + randPoint(params.R0,1,'extreme');
        else
            params.x0 = x0_m + randPoint(params.R0,1);
        end
        if rand(1) < 0.3
            params.u = u_m + randPoint(params.U, N_k,'extreme');
        else
            params.u = u_m + randPoint(params.U, N_k);
        end
        % simulate the system dynamics to obtain the output trajectory
        [~,~,~,y_m{s}] = simulate(sys,params);
    end
    % create testCase for each nominal trajectory
    testSuite{m} = testCase(y_m, u_m', x0_m, sys.dt);
end
end


% ------------------------------ END OF CODE ------------------------------
