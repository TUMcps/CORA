 function res = test_nonlinear_linError()
% test_nonlinear_linError - test if the linearization error for nonlinear
%                           systems is computed correctly
%
% Syntax:  
%    res = test_nonlinear_linError()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      12-November-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

    res = 1;

    % set options for reachability analysis
    options.tensorOrder = 2;
    options.taylorTerms = 5;

    options.tStart = 0;
    options.tFinal = 0.0021;
    options.timeStep = 0.0021;                  

    options.zonotopeOrder = 50;
    options.reductionTechnique = 'girard';

    options.maxError = 1e20 * ones(18,1);
    options.reductionInterval = Inf;

    options.U = zonotope([0;0]);
    options.uTrans = [0;0];

    options.advancedLinErrorComp = 0;

    x0 = zeros(18,1);
    x0(4) = 19.9536;
    x0(5) = 18.6195;
    x0(8) = 0.7098;

    G =   [[0 0 0 0.0203 0.0125 0.0123 0 0 0 0; ...
            0 0 0 -0.0068 -0.0015 0.0001 0.0071 0 0 0; ...
            0 0 0 0.1635 0.1524 -0.0041 0 0 0 0; ...
            0.15 0.1443 0.1391 0 0 0 0 0 0 0.06; ...
            -0.0016 -0.0045 -0.0072 0 0 0 0 0 0.0796 -0.0037; ...
            0 0 0 -0.0103 -0.0022 -0.0118 -0.0089 0.0786 0 0; ...
            0 0 0 0.0006 -0.0634 0 0 0 0 0; ...
            -7.1991 -0.0594 -0.0117 0 0 0 0 0 0 -0.0099]; eye(10)];

    options.R0 = zonotope([x0,G]);


    % construct the nonlinear system
    sys = nonlinearSys(18,2,@highorderBicycleDynamics,options);

    % compute the reachable set
    [~,R] = reach(sys,options);
    
    % extract the set of linearization errors
    err = R{1}{1}.error;
    linError = zonotope([zeros(length(err),1),diag(err)]);
    
    % evaluate the linearization error for a set of random points
    p.u = center(options.U) + options.uTrans;
    f0prev = sys.mFile(0,center(options.R0),p.u);
    p.x = center(options.R0) + f0prev*0.5*options.timeStep;

    f0 = sys.mFile(0,p.x,p.u);
    [A,~] = sys.jacobian(p.x,p.u);
    
    N = 10000;
    points = zeros(length(x0),N);

    for i = 1:N
        p = randPoint(options.R0);

        temp = sys.mFile(0,p,[0;0]);
        points(:,i) = temp - A * (p - sys.linError.p.x) - f0;
    end
    
%     % visualize the result
%     plot(points(1,:),points(3,:),'.k');
%     hold on
%     plot(linError,[1,3],'r');
    
    % check if the set of linerization error contains all randomly computed
    % points
    linError = interval(linError);
    
    linError = linError([1 3 5 6]);
    points = points([1 3 5 6],:);
    
    for i = 1:N
       if ~in(linError,points(:,i))
          error('test_nonlinear_error failed!'); 
       end
    end

%------------- END OF CODE --------------