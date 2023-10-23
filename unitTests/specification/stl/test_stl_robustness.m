function res = test_stl_robustness()
% test_stl_robustness - unit test function of robustness
%
% Syntax:
%    res = test_stl_robustness
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       03-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    res = true;

    % Analytical Tests ----------------------------------------------------

    % analytical tests for the until-operator

    % trajectory
    phi = -pi/2:0.01:0;
    x = [cos(phi'),sin(phi')];
    t = linspace(0,1,length(phi))';
    
    % STL formula
    y = stl('y',2);
    point = sqrt(2)/2 - 0.01;
    eq = until(y(2) < -point,y(1) > point,interval(0,1));
    
    % model checking
    if robustness(eq,x,t) < 0
        throw(CORAerror('CORA:testFailed'));
    end

    % slightly change STL formula
    y = stl('y',2);
    point = sqrt(2)/2 + 0.01;
    eq = until(y(2) < -point,y(1) > point,interval(0,1));

    % model checking
    if robustness(eq,x,t) > 0
        throw(CORAerror('CORA:testFailed'));
    end


    % analytical tests for the globally-operator

    % trajectory
    x = linspace(0,1,100)';
    t = linspace(0,1,100)';

    % STL formula
    y = stl('y',1); 
    eq = globally(y(1) > 0.1,interval(0,1));

    % model checking
    if robustness(eq,x,t) > 0
        throw(CORAerror('CORA:testFailed'));
    end

    % different STL formula
    eq = globally(y(1) >= 0.2 & y(1) <= 0.3, interval(0.2,0.3));

    % model checking
    if robustness(eq,x,t) < 0
        throw(CORAerror('CORA:testFailed'));
    end   


    % analytical tests for the finally-operator

    % trajectory
    x = linspace(0,1,100)';
    t = linspace(0,1,100)';

    % STL formula
    y = stl('y',1); 
    eq = finally(y(1) > 1.1,interval(0,1));

    % model checking
    if robustness(eq,x,t) > 0
        throw(CORAerror('CORA:testFailed'));
    end

    % different STL formula
    eq = finally(y(1) >= 0.25, interval(0.2,0.3));

    % model checking
    if robustness(eq,x,t) < 0
        throw(CORAerror('CORA:testFailed'));
    end   


    % analytical tests for different set representations

    % dynamical system
    sys = linearSys([0 -1;1 0],[0;0]);

    % compute reachable set
    params.R0 = zonotope([0;-1] + 0.1*interval([-1;-1],[1;1]));
    params.tFinal = 2;

    options.timeStep = 0.1;
    options.taylorTerms = 5;
    options.zonotopeOrder = 10;

    R = reach(sys,params,options);

    % simulate single trace
    simOpts.x0 = center(params.R0);
    simOpts.tFinal = params.tFinal;

    [t_,x_] = simulate(sys,simOpts);

    t = 0:0.01:params.tFinal;
    [~,ind] = unique(t_);
    x = [interp1(t_(ind),x_(ind,1),t,'linear','extrap'); ...
         interp1(t_(ind),x_(ind,2),t,'linear','extrap')]';

    % test for halfspace predicates
    y = stl('y',2);
    eq = until(y(2) < -0.697,y(1) > 0.697,interval(0,2));

    robReach = robustness(eq,R);
    robSim = robustness(eq,x,t);

    if robSim < robReach
        throw(CORAerror('CORA:testFailed'));
    end

    % test for polytope predicates
    y = stl('y',2);
    eq = globally(y(1) + 0.1*y(2) > -0.5 | ...
                        y(2) - 0.05*y(1) < 0.5, interval(0,2));

    robReach = robustness(eq,R);
    robSim = robustness(eq,x,t);

    if robSim < robReach
        throw(CORAerror('CORA:testFailed'));
    end

    % test for level set predicates
    y = stl('y',2);
    eq = globally(sqrt(y(1)^2 + y(2)^2) > 0.8 | ...
                        sin(y(2)) - 0.05*y(1) < 0.5, interval(0,2));

    robReach = robustness(eq,R);
    robSim = robustness(eq,x,t);

    if robSim < robReach
        throw(CORAerror('CORA:testFailed'));
    end


    % Random Tests --------------------------------------------------------

    % check if robustness is consistent with model checking results
    for i = 1:5

        % generate random linear system
        sys = linearSys(-1+2*rand(2,2),-1+2*rand(2,1));
    
        % generate random trace
        simOpts.x0 = -2+4*rand(2,1);
        simOpts.u = -2 + 4*rand(1,10);
        simOpts.tFinal = 2;
    
        [t_,x_] = simulate(sys,simOpts);
    
        t = 0:0.02:2;
        [~,ind] = unique(t_);
        x = [interp1(t_(ind),x_(ind,1),t,'linear','extrap'); ...
             interp1(t_(ind),x_(ind,2),t,'linear','extrap')]';
        
        % generate random STL formula
        dom = interval(min(x,[],1)',max(x,[],1)');
        eq = stl.generateRandom('Dimension',2,'FinalTime',1, ...
                                            'TimeStep',0.1,'Domain',dom);
    
        % model check formula
        resOrig = modelCheckTrace(negationNormalForm(eq),x,t);

        % compute robustness
        rob = robustness(eq,x,t);

        % check if results are consistent
        if ((resOrig && rob < 0) || (~resOrig && rob > 0)) && ~isinf(rob)
            throw(CORAerror('CORA:testFailed'));
        end
    end

    % check if robustness for reachable set is consistent with robustness
    % for traces contained inside the reachable set
    for i = 1:5

        % generate random linear system
        sys = linearSys(-1+2*rand(2,2),-1+2*rand(2,1));
    
        % generate reachability parameter
        params.R0 = zonotope((-2+4*rand(2,1)) + interval([-1;-1],[1;1]));
        params.u = -2 + 4*rand(1,20);
        params.tFinal = 1;

        % compute reachable set
        options.timeStep = 0.05;
        options.taylorTerms = 5;
        options.zonotopeOrder = 10;

        R = reach(sys,params,options);

        % simulate random trajectories inside the reachable set
        xsim = {};

        for j = 1:10
            simOpts.x0 = randPoint(params.R0);
            simOpts.u = params.u;
            simOpts.tFinal = params.tFinal;
    
            [t_,x_] = simulate(sys,simOpts);
    
            t = 0:0.01:1;
            [~,ind] = unique(t_);
            x = [interp1(t_(ind),x_(ind,1),t,'linear','extrap'); ...
                 interp1(t_(ind),x_(ind,2),t,'linear','extrap')]';
            xsim{end+1} = x;
        end
        
        % generate random STL formula
        dom = [];

        for j = 1:length(R.timeInterval.set)
            dom = dom | interval(R.timeInterval.set{j});
        end

        eq = stl.generateRandom('Dimension',2,'FinalTime',1, ...
                                            'TimeStep',0.1,'Domain',dom);

        % compute robustness for the reachable set
        robReach = robustness(eq,R);

        % compute robustness for the simulated traces
        robSim = inf;

        for j = 1:length(robSim)
            robSim = min(robSim,robustness(eq,xsim{j},t));
        end

        % check if results are consistent
        if robSim < robReach
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
