function sys = identify(traj, options)
% identify - Identify a nonlinear ARX system from trajectory data using
%   genetic programming
%
% Syntax:
%    sys = nonlinearARX.identify(traj,options)
%
% Inputs:
%    traj - object of class "trajectory" storing the trajectory data
%    options - algorithm options for system identification
%
% Outputs:
%    sys - identified nonlinearARX object
%
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-Conformant System
%        Identification," arXiv, 2025. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/identify

% Authors:       Laura Luetzow
% Written:       03-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse the input arguments 
params = options.params;
options = rmfield(options,'params');

% create nonlinear ARX object and validate the input specifications
dim_y = size(traj(1).y, 1);
dim_u = size(traj(1).u, 1);
dt = traj(1).dt;
n_p = options.id.p;
dummySys = nonlinearARX([],dt,dim_y,dim_u,n_p);
[params,options] = validateOptions(dummySys,params,options);

% reformat training and validation data
[xtrain,ytrain] = aux_testSuite2regress(traj, options.id.p, dt);
[xval,yval] = aux_testSuite2regress(options.id.testSuite_val, options.id.p, dt);

if options.id.verbose
    fprintf("Approximate dynamics with %s. \n", options.alg);
end

% call the selected algorithm
switch options.idAlg
    case 'gp'
        f = aux_identifyGP(xtrain, ytrain, xval, yval, params, options);
    case 'cgp'
        f = aux_identifyCGP(xtrain, ytrain, xval, yval, params, options);
end

% transform nonlinear function to sys object
sys = nonlinearARX(f,dt, dim_y,dim_u, n_p);

end


% Auxiliary functions -----------------------------------------------------

function [f, fitness_sum] = aux_identifyGP(xtrain, ytrain, xval, yval, params, options)
% normal genetic programming

options.id.xval = xval;
options.id.xtrain = xtrain;
fitness_sum = 0;
func = "@(y,u) [";
for i_y = 1:size(ytrain,2)
    options.id.ytrain = ytrain(:,i_y);
    options.id.yval = yval(:,i_y);
    file_approx_iy = options.id.filename + sprintf("_dim%d",i_y);

    % start new gp run
    if options.id.verbose
        fprintf('Output dimension %d. \n', i_y)
    end
    timerVal = tic;
    gp = rungp(@(x) config_gp(x, params, options, 'gp'));
    T = toc(timerVal);
    if options.id.save_res
        save(file_approx_iy,"gp", "T");
    end
    expr = gppretty(gp,'valbest');
    func = func + string(expr) + ";";
    fitness_sum = fitness_sum + gp.results.best.fitness;
end
f = eval(func + "]");
end


function f = aux_identifyCGP(xtrain, ytrain, xval, yval, params, options)
% genetic programming with conformance cost

options.id.xval = xval;
options.id.xtrain = xtrain;
options.id.yval = yval;
options.id.ytrain = ytrain;
options.id.pop_pre = true;
timerVal = tic;
% create initial population with normal genetic programming
del_res = false;

if ~isfield(options.id, 'cgp_file_pop_pre') && ...
        options.id.gp_num_gen > options.id.cgp_num_gen
    % create initial propulation with normal genetic programming
    options_gp = options;
    options_gp.id.gp_num_gen = options.id.gp_num_gen - options.id.cgp_num_gen;
    options_gp.id.gp_runs = 1;
    options_gp.id.save_res = true;
    [~, fitness] = aux_identifyGP(xtrain, ytrain, xval, yval, params, options_gp);
    options.id.cgp_file_pop_pre = options.id.filename;
    options.id.cgp_conf_value = 5*fitness;
    if ~options.id.save_res
        del_res = true;
    end
end

% use validation testSuite
if ~isfield(options.id, 'testSuite_val')
    params.testSuite = options.id.testSuite_val;
end
gp = rungp(@(x)config_gp(x, params, options, 'cgp'));
T = toc(timerVal);
if options.id.save_res
    save(options.id.filename,"gp", "T");
end
if del_res
    % delete intermediary results
    for i_y = 1:size(ytrain,2)
        file_approx_iy = options.id.filename + sprintf("_dim%d",i_y);
        delete(file_approx_iy+".mat");
    end
end

expr = gppretty(gp,'best');
func = "@(y,u) [";
for i_y = 1:size(ytrain,2)
    func = func + string(expr(i_y)) + ";";
end
f = eval(func + "]");
end


function [x,y] = aux_testSuite2regress(testSuite, p, dt)
% x = [y_1(k-p) y_2(k-p) ... y_n(k-p) y_1(k-1) y_2(k-1) ... y_n(k-1) ...
%      u_1(k-p) u_2(k-p) ... u_n(k-p) u_1(k) u_2(k) ... u_n(k)]

if ~testSuite(1).constant_dt
    % convert data to uniform time step
    testSuite = uniformTimeStepSize(testSuite, dt);
end

total_size = length(testSuite) * (size(testSuite(1).y,2)-p) * size(testSuite(1).y, 3);
x = zeros(total_size, size(testSuite(1).y, 1)*p + size(testSuite(1).u, 1)*(p+1));
y = zeros(total_size, size(testSuite(1).y, 1));
index = 1;
for m = 1:length(testSuite)
    y_m = testSuite(m).y;
    u_m = testSuite(m).u;
    for k = p+1:size(y_m,2)
        x_k = [reshape(y_m(:,k-p:k-1,:), 1, [], size(y_m, 3)) ...
            repmat(reshape(u_m(:,k-p:k), 1, [], 1), 1, 1, size(y_m, 3))];
        y_k = y_m(:,k,:);
        x(index:index+size(y_k, 3)-1, :) = squeeze(x_k)';
        y(index:index+size(y_k, 3)-1, :) = squeeze(y_k)';
        index = index + size(y_k, 3);
    end
end
end


% ------------------------------ END OF CODE ------------------------------
