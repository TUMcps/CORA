function res = testLong_nonlinearSysDT_tensorCreation()
% testLong_nonlinearSysDT_tensorCreation - unit test function
%    for the creation of the third-order-tensor file
%
% Checks different scenarios of settings, where each scenario results in a
% different third-order tensor
%
% Syntax:
%    res = testLong_nonlinearSysDT_tensorCreation()
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

% Options -----------------------------------------------------------------

% time
params.tFinal = 0.03; 

% initial set
params.R0 = zonotope([[-0.15;-45],diag([0.005;3])]); % initial set

% algorithm settings
options.zonotopeOrder=50; % maximum zonotope order
options.tensorOrder = 3;
options.errorOrder = 10;

% system inputs
params.U = zonotope([zeros(2,1),diag([0.1;2])]);


% Test Cases --------------------------------------------------------------

% For each of the considered scenarios, a different third-order-tensor file
% is created

for i = 1:2
    for j = 1:2
        for h = 1:2
            for m = 1:2

                options_ = options;

                % replacements
                if i == 1
                    options_.lagrangeRem.replacements = @(x,u) exp(-8750/(x(2) + 350)); 
                end

                % parallel execution
                if j == 1
                    options_.lagrangeRem.tensorParallel = true;
                end

                % taylor models 
                if h == 1
                    options_.lagrangeRem.method = 'taylorModel';
                end
                
                % tensor order
                if m == 1
                   options_.tensorOrder = 2; 
                else
                   options_.tensorOrder = 3; 
                end

                % create system 
                dt = 0.015;
                fun = @(x,u) cstrDiscr(x,u,dt);
                sys = nonlinearSysDT('stirredTankReactor',fun,dt);

                % compute reachable set
                reach(sys, params, options_);
            end
        end           
    end
end

% test is successfull if no error occured during execution
res = true;

% ------------------------------ END OF CODE ------------------------------
