function res = test_nonlinearDT_tensorCreation()
% test_nonlinearDT_tensorCreation - unit_test_function for the creation of
%                                   the third-order-tensor file
%
% Checks different scenarios of settings, where each scenario results in a
% different third-order tensor
%
% Syntax:  
%    res = test_nonlinearDT_tensorCreation()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
% 
% Author:       Niklas Kochdumper
% Written:      02-August-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

res = 0;
dim=6;


% Options -----------------------------------------------------------------

% time
options.tStart = 0; % start time
options.tFinal = 0.02; % final time
options.timeStep = 0.01; %time step size 

% initial set
options.x0 = [-0.15;-45]; % center of initial set 
options.R0 = zonotope([options.x0,diag([0.005;3])]); % initial set

% algorithm settings
options.taylorTerms=4; % number of taylor terms 
options.zonotopeOrder=50; % maximum zonotope order
options.tensorOrder = 3;
options.errorOrder = 10;
options.reductionTechnique='girard';

% system inputs
options.uTrans = 0;
options.U = zonotope([zeros(2,1),diag([0.1;2])]);


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
                    options_.replacements = @(x,u) exp(-8750/(x(2) + 350)); 
                end

                % parallel execution
                if j == 1
                    options_.tensorParallel = 1;
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
                sys = nonlinearSysDT(2,2,@cstrDiscr,options_);

                % compute reachable set
                reach(sys, options_);
            end
        end           
    end
end

% test is successfull if no error occured during execution
res = 1;

%------------- END OF CODE --------------