function res = test_nonlinearParam_tensorCreation()
% test_nonlinearParam_tensorCreation - unit_test_function for the creation
%                                      of the third-order-tensor file
%
% Checks different scenarios of settings, where each scenario results in a
% different third-order tensor
%
% Syntax:  
%    res = test_nonlinearParam_tensorCreation()
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
options.tStart=0; % start time
options.tFinal=8; % final time
options.timeStep=4; %time step size 

% initial set
options.x0=[2; 4; 4; 2; 10; 4]; % center of initial set 
options.R0=zonotope([options.x0,0.2*eye(dim)]); % initial set

% algorithm settings
options.taylorTerms=4; % number of taylor terms 
options.zonotopeOrder=50; % maximum zonotope order
options.intermediateOrder=5;
options.originContained = 0;
options.advancedLinErrorComp = 1;
options.tensorOrder = 3;
options.reductionTechnique='girard';
options.errorOrder=1;
options.reductionInterval=1e3;
options.maxError = 1*ones(dim,1);

% parameter values
options.paramInt=0.015;

% system inputs
options.uTrans = 0;
options.U = zonotope([0,0.005]); %input for reachability analysis


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
                        options_.replacements = @(x,u,p) (1.661042594276257564748421423659*p(1))/x(1)^(5/2);
                    end

                    % parallel execution
                    if j == 1
                        options_.tensorParallel = 1;
                    end

                    % reachability algorithm
                    if k == 1
                        options_.advancedLinErrorComp = 1;
                    else
                        options_.advancedLinErrorComp = 0;
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
                    sys = nonlinParamSys(6,1,1,@tank6paramEq,options_);

                    % compute reachable set
                    reach(sys, options_);
                end
            end           
        end
    end
end

% test is successfull if no error occured during execution
res = 1;