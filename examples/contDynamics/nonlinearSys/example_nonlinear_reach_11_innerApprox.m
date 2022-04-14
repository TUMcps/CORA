function example_nonlinear_reach_11_innerApprox()
% example_nonlinear_reach_11_innerApprox - example demonstrating the 
%       computation of inner-approximation of reachable sets with the
%       algorithms from [1] and [2]. The example can be found in Sec. 4.1 
%       of [2].
%
%
% Syntax:  
%    example_nonlinear_reach_11_innerApprox
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% References:
%    [1] N. Kochdumper and M. Althoff. "Computing Non-Convex Inner-
%        Approximations of Reachable Sets for Nonlinear Continuous Systems"
%        CDC 2020
%    [2] E. Goubault and S. Putot. "Forward Inner-Approximated Reachability
%        of Non-Linear Continuous Systems", HSCC 2017 
%    [3] E. Goubault and S. Putot. "Robust Under-Approximations and 
%        Application to Reachability of Non-Linear Control Systems With 
%        Disturbances", Control System Letters 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachInner

% Author:       Niklas Kochdumper
% Written:      21-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


    % Parameter -----------------------------------------------------------

    params.tFinal = 1;
    params.R0 = interval([0.9;0],[1;0.1]);


    % Reachability Settings -----------------------------------------------

    % settings for the algorithm from [1]
    options1.algInner = 'scale';
    
    options1.splits = 2;
    options1.iter = 2;
    options1.orderInner = 5;
    options1.scaleFac = 0.95;
    
    options1.timeStep = 0.001;                           
    options1.taylorTerms = 10;                            
    options1.zonotopeOrder = 50;       
    options1.intermediateOrder = 20;
    options1.errorOrder = 10;

    % settings for the algorithm from [2]
    options2.algInner = 'proj';
    options2.timeStep = 0.05;
    options2.taylorOrder = 3;
    options2.taylmOrder = 10;
    
    % settings for the algorithm from [3]
    options3.algInner = 'parallelo';
    options3.alg = 'lin';
    options3.timeStep = 0.05;
    options3.zonotopeOrder = 20;
    options3.taylorTerms = 10;
    options3.tensorOrder = 2;


    % System Dynamics -----------------------------------------------------

    brusselator = @(x,u) [1-2*x(1) + 3/2 * x(1)^2*x(2); ...
                          x(1)-3/2*x(1)^2*x(2)];
    
    sys = nonlinearSys(brusselator); 


    % Reachability Analysis -----------------------------------------------
    
    % compute inner-approximation with the algorithm from [1] 
    tic
    [Rin1,Rout1] = reachInner(sys,params,options1);
    tComp = toc;
    
    disp(['Computation time (scale): ',num2str(tComp),' s']);

    % compute inner-approximation with the algorithm from [2] 
    tic
    [Rin2,Rout2] = reachInner(sys,params,options2);
    tComp = toc;
    
    disp(['Computation time (proj): ',num2str(tComp),' s']);
    
    % compute inner-approximation with the algorithm from [3] 
    tic
    [Rin3,~] = reachInner(sys,params,options3);
    tComp = toc;
    
    disp(['Computation time (parallelo): ',num2str(tComp),' s']);
    
    
    % Visualization -------------------------------------------------------

    % compare inner with outer approximation
    figure; hold on;
    hanOut = plot(Rout1.timePoint.set{end},[1,2],'r');
    hanInScal = plot(Rin1.timePoint.set{end},[1,2],'g');
    hanInProj = plot(Rin2.timePoint.set{end},[1,2],'b');
    hanInPar = plot(Rin3.timePoint.set{end},[1,2],'c');
    
    l = legend([hanOut,hanInScal,hanInProj,hanInPar],'outer approx.', ...
               'inner approx. (scal)','inner approx. (proj)', ...
               'inner approx. (parallelo)');
    set(l,'Location','northwest');
    xlabel('x_1');
    ylabel('x_2');
    
    % plot inner approximation over time
    figure; hold on
    plotOverTime(Rout2,1,'FaceColor',[.7 .7 .7],'EdgeColor','none');
    hOut = plotOverTime(Rout2,2,'FaceColor',[.7 .7 .7],'EdgeColor','none');
    hIn1 = plotOverTime(Rin2,1,'FaceColor','none','EdgeColor','r');
    hIn2 = plotOverTime(Rin2,2,'FaceColor','none','EdgeColor','b');
    
    xlabel('time');
    ylabel('x_1, x_2')
    box on
    l = legend([hOut;hIn1;hIn2],'Outer approx.','Inner approx. (x_1)', ...
                                'Inner approx. (x_2)');
    set(l,'Location','northeast');
end

%------------- END OF CODE -------------