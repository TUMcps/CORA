function res = test_nonlinearSys_simplifyOptimize
% test_nonlinearSys_simplifyOptimize - unit test function for the creation
%    of tensors with setting options.simplify = 'optimize'; Checks if the
%    results with options.simplify = 'optimize' are identical to the
%    results without optimization
%
% Syntax:
%    res = test_nonlinearSys_simplifyOptimize
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       12-March-2020
% Last update:   10-October-2024 (MW, use derivatives instead of reach)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% general parameters
N = 10;
R = interval([-925;-425;-1;-1;0],[-875;-375;1;1;10]);

% Hessian Tensors ------------------------------------------------------

options.tensorOrder = 2;
options.tensorOrderOutput = 2;
options.lagrangeRem.tensorParallel = false;
options.lagrangeRem.method = 'interval';

% create nonlinear system objects
sys1 = nonlinearSys('testfun1',@aux_testDyn);
sys2 = nonlinearSys('testfun2',@aux_testDyn);

% call derivatives to create tensors
options.lagrangeRem.simplify = 'none';
derivatives(sys1,options);
setHessian(sys1,'int');

options.lagrangeRem.simplify = 'optimize';
derivatives(sys2,options);
setHessian(sys2,'int');

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
        assertLoop(isequal(H1{j},H2{j},1e-14),i,j)
    end
end


% Third-Order Tensors -----------------------------------------------------

% create nonlinear system objects anew
sys1 = nonlinearSys('testfun1',@aux_testDyn);
sys2 = nonlinearSys('testfun2',@aux_testDyn);

% set tensor order to 3 to create third-order tensor
options.tensorOrder = 3;

% call derivatives to create tensors
options.lagrangeRem.simplify = 'none';
derivatives(sys1,options);
setThirdOrderTensor(sys1,'int');

options.lagrangeRem.simplify = 'optimize';
derivatives(sys2,options);
setThirdOrderTensor(sys2,'int');

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
            assertLoop(isequal(T1{j,k},T2{j,k},1e-14),i,j,k)
        end
    end
end

% test completed
res = true;
    
end


% Auxiliary functions -----------------------------------------------------

function dx = aux_testDyn(x,u)

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
