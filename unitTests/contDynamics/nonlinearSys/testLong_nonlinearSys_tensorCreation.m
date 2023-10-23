function res = testLong_nonlinearSys_tensorCreation
% testLong_nonlinearSys_tensorCreation - unit_test_function for the
%    creation of the third-order-tensor file
%
%    Checks different scenarios of settings, where each scenario results in
%    a different third-order tensor
%
% Syntax:
%    res = testLong_nonlinearSys_tensorCreation()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       02-August-2018
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;
dim_x = 6;

% Parameters --------------------------------------------------------------

params.tFinal=8; % final time
params.R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]); % initial set
params.U = zonotope([0,0.005]); %input for reachability analysis


% Reachability Settings ---------------------------------------------------

options.timeStep=4; %time step size 
options.taylorTerms=4; % number of taylor terms 
options.zonotopeOrder=50; % maximum zonotope order
options.intermediateOrder=5;
options.alg = 'poly';
options.tensorOrder = 3;
options.errorOrder=1;
options.reductionInterval=1e3;
options.maxError = 1*ones(dim_x,1);


% Test Cases --------------------------------------------------------------

% For each of the considered scenarios, a different third-order-tensor file
% is created

for i = 1:2
    for j = 1:2
        for k = 1:2
            for h = 1:2
                for m = 1:2
            
                    options_ = options;

                    % replacements
                    if i == 1
                        options_.lagrangeRem.replacements = @(x,u) 897680497035489/(36028797018963968*x(5)^(5/2)); 
                    end

                    % parallel execution
                    if j == 1
                        options_.lagrangeRem.tensorParallel = true;
                    end

                    % reachability algorithm
                    if k == 1
                        options_.alg = 'poly';
                    else
                        options_.alg = 'lin';
                    end

                    % taylor models 
                    if h == 1
                        options_.lagrangeRem.method = 'taylorModel';
                    end
                    
                    % tensor order
                    if m == 1 && ~strcmp(options.alg,'poly')
                        options_.tensorOrder = 2;
                    else
                        options_.tensorOrder = 3;
                    end

                    % create system 
                    sys = nonlinearSys(@tank6Eq);

                    % compute reachable set
                    reach(sys, params, options_);
                end
            end           
        end
    end
end

% test is successful if no error occured during execution
res = true;

% ------------------------------ END OF CODE ------------------------------
