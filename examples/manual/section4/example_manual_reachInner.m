function example_manual_reachInner()
% example_manual_reachInner - example from the manual demonstrating the 
%                             inner-approximation of reachable sets
%
% Syntax:  
%   example_manual_reachInner()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_manual_reach
%
% Author:        Niklas Kochdumper
% Written:       18-May-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % system dynamics
    f = @(x,u) [1-2*x(1) + 3/2 * x(1)^2*x(2); ...
                x(1)-3/2*x(1)^2*x(2)];

    sys = nonlinearSys(f); 

    % parameter
    params.tFinal = 1;
    params.R0 = interval([0.75;0],[1;0.25]);

    % reachability settings
    options.algInner = 'scale';
    options.timeStep = 0.001;                           
    options.taylorTerms = 10;                            
    options.zonotopeOrder = 50;       
    options.intermediateOrder = 20;
    options.errorOrder = 10;

    % reachability analysis
    [Rin,Rout] = reachInner(sys,params,options);

    % visualization
    figure; hold on; box on;
    plot(Rout.timePoint.set{end},[1,2],'b');
    plot(Rin.timePoint.set{end},[1,2],'r');
    axis equal
    xlim([0.6,0.8])
    ylim([0.45,0.65])
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    
%------------- END OF CODE --------------