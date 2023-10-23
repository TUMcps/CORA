function res = testLong_nonlinParamSys_tensorCreation
% testLong_nonlinParamSys_tensorCreation - unit test function for
%    the creation of the third-order-tensor file
%
% Checks different scenarios of settings, where each scenario results in a
% different third-order tensor
%
% Syntax:
%    res = testLong_nonlinParamSys_tensorCreation
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

% Parameters --------------------------------------------------------------

params.R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(6)]);
params.U = zonotope([0,0.005]); 
params.tFinal = 8; 
params.paramInt = 0.015;


% Reachability Settings ---------------------------------------------------

options.timeStep = 4; 
options.taylorTerms = 5;
options.zonotopeOrder = 50;
options.intermediateTerms = 5;
options.errorOrder = 1;
options.intermediateOrder = 1;

options.alg = 'poly';
options.tensorOrder = 3;


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
                        options_.lagrangeRem.replacements = @(x,u,p) (1.661042594276257564748421423659*p(1))/x(1)^(5/2);
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
                    if m == 1 && ~strcmp(options_.alg,'poly')
                        options_.tensorOrder = 2;
                    else
                        options_.tensorOrder = 3; 
                    end

                    % create system 
                    sys = nonlinParamSys(@tank6paramEq);

                    % compute reachable set
                    reach(sys, params, options_);
                end
            end           
        end
    end
end

% test is successful if no error occurred during execution
res = true;

% ------------------------------ END OF CODE ------------------------------
