function res = testLong_nonlinearSys_simplifyOptimize
% testLong_nonlinearSys_simplifyOptimize - unit test function for the creation
%    of tensors with setting options.simplify = 'optimize'
%
% Checks if the results with options.simplify = 'optimize' are identical to
% the results without optimization
%
% Syntax:
%    res = testLong_nonlinearSys_simplifyOptimize
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       12-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % general parameter
    N = 100;
    R = interval([-925;-425;-1;-1;0],[-875;-375;1;1;10]);
    
    % parameter for reachability analysis
    params.tFinal = 0.01;
    params.R0 = zonotope(R);
    
    % settings for reachability analysis
    options.alg = 'lin';
    options.zonotopeOrder = 10;
    options.taylorTerms = 5;
    options.timeStep = 0.01;
    

    % Hessian Tensor ------------------------------------------------------
    
    options_ = options;
    options_.tensorOrder = 2;
    
    % create nonlinear system objects
    sys1 = nonlinearSys(@aux_testDyn1);
    sys2 = nonlinearSys(@aux_testDyn2);
    
    % call reach so that tensors are created
    reach(sys1,params,options_);
    
    options_.lagrangeRem.simplify = 'optimize';
    reach(sys2,params,options_);
    
    % check for multiple random points if the resulting tensor is identical
    for i = 1:N
       
        % draw random point
        x = interval(randPoint(R));
        u = interval(1);
        
        % evaluate tensors
        H1 = sys1.hessian(x,u);
        H2 = sys2.hessian(x,u);
        
        % compare the results
        for j = 1:length(H1)
           
            temp1 = center(H1{j});
            temp2 = center(H2{j});
            
            if ~all(all(withinTol(temp1,temp2,1e-14)))
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
    
    
    % Third-Order Tensor --------------------------------------------------

    options_ = options;
    options_.errorOrder = 10;
    options_.intermediateOrder = 10;
    options_.tensorOrder = 3;
    
    % create nonlinear system objects
    sys1 = nonlinearSys(@aux_testDyn1);
    sys2 = nonlinearSys(@aux_testDyn2);
    
    % call reach so that tensors are created
    reach(sys1,params,options_);
    
    options_.lagrangeRem.simplify = 'optimize';
    reach(sys2,params,options_);
    
    % check for multiple random points if the resulting tensor is identical
    for i = 1:N
       
        % draw random point
        x = interval(randPoint(R));
        u = interval(1);
        
        % evaluate tensors
        T1 = sys1.thirdOrderTensor(x,u);
        T2 = sys2.thirdOrderTensor(x,u);
        
        % compare the results
        for j = 1:size(T1,1)
            for k = 1:size(T1,2)
           
                temp1 = center(T1{j,k});
                temp2 = center(T2{j,k});

                if ~all(all(withinTol(temp1,temp2,1e-14)))
                    throw(CORAerror('CORA:testFailed'));
                end
            end
        end
    end
    
    res = true;
    
end


% Auxiliary functions -----------------------------------------------------

function dx = aux_testDyn1(x,u)

    dx(1,1) = x(3);
    dx(2,1) = x(4);
    dx(3,1) = 0.0002624*x(2) - 0.576038456806855922248709263104*x(1) - ...
              19.2299796*x(3) + 0.008750586984672*x(4) - ...
              (1434960000000000000.0*(x(1) + 42164000))/((x(1) + ...
              42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;
    dx(4,1) = - 0.0002624*x(1) - 0.575980856806855922248709263104*x(2) - ...
              0.008750586984672*x(3) - 19.2299766*x(4) - ...
              (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);
    dx(5,1) = 1;
end

function dx = aux_testDyn2(x,u)

    dx(1,1) = x(3);
    dx(2,1) = x(4);
    dx(3,1) = 0.0002624*x(2) - 0.576038456806855922248709263104*x(1) - ...
              19.2299796*x(3) + 0.008750586984672*x(4) - ...
              (1434960000000000000.0*(x(1) + 42164000))/((x(1) + ...
              42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;
    dx(4,1) = - 0.0002624*x(1) - 0.575980856806855922248709263104*x(2) - ...
              0.008750586984672*x(3) - 19.2299766*x(4) - ...
              (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);
    dx(5,1) = 1;
end

% ------------------------------ END OF CODE ------------------------------
