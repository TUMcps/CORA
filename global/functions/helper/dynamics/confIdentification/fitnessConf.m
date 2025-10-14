function fitness = fitnessConf(gp, id_curr)
% fitnessConf - compute the fitness of a given gp function using
%                       conformance identification
%
% Syntax:
%    fitness = fitnessConf(gp, id_curr)
%
% Inputs:
%    gp - genetic programming object
%    id_curr - identifier for the individual, which is to be evaluated
%
% Outputs:
%    fitness - fitness value
%
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System
%        Identification," arXiv, 2024. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Laura Luetzow
% Written:       08-March-2024             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 2
    id_curr = gp.state.current_individual;
end

% load data
params_id_init = gp.userdata.params_conf;
options = gp.userdata.options_conf;
n_m = min([gp.userdata.n_m_conf, length(params_id_init.testSuite)]); 
n_id = floor(length(params_id_init.testSuite)/n_m);

% create dynamical system object
try
    expr = gppretty(gp,id_curr);
    func = "@(y,u) [";
    for i_y = 1:size(gp.userdata.ytrain,2)
        func = func + string(expr(i_y)) + ";"; %runtree(gp,'valbest');
    end
    f = eval(func + "]");
    sys_approx = nonlinearARX("NARX_gpcTrain_"+string(id_curr),f,...
        gp.userdata.dt,size(params_id_init.testSuite(1).y,1), ...
        size(params_id_init.testSuite(1).u,1), gp.userdata.p);

    % compute conformance cost
    timerVal = tic;
    fval = 0;
    meas_dist = zeros(n_id,1);

    testSuite = params_id_init.testSuite;
    for i=0:n_id-1
        params_id_init.testSuite = testSuite(n_m*i+1:n_m*(i+1));
        for j=1:n_m
            % maximum distance between measurements of the same test case
            meas_dist(i+1) = meas_dist(i+1) + sum(max(abs(diff(params_id_init.testSuite(j).y,3)),[],3),'all');
        end
        [~, results] = conform(sys_approx,params_id_init,options);
        fval = fval + results.fval/meas_dist(i+1);
    end
catch ME
    display(ME.message)
    fval = inf;
end
T=toc(timerVal);
fitness = mean(meas_dist)*fval + round(T,1);
end

% ------------------------------ END OF CODE ------------------------------
