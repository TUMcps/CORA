function res = test_stl_modelCheckReachSet()
% test_stl_modelCheckReachSet - unit test function of modelCheckReachSet
%
% Syntax:
%    res = test_stl_modelCheckReachSet
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
% Written:       09-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = true;
    alg = {'rtl','sampledTime','signals'};

    % Analytical Test 1 ---------------------------------------------------

    % analytical test for the until-operator

    % dynamic system, parameters, and options
    sys = linearSys([0 -1; 1 0],[0;0]);

    params.R0 = zonotope([0;-1]);
    params.tFinal = 2;

    options.timeStep = 0.5;
    options.zonotopeOrder = 10;
    options.taylorTerms = 10;

    % STL formula
    y = stl('y',2);
    point = sqrt(2)/2;
    eq = until(y(2) < -point,y(1) > point,interval(0,2));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end


    % Analytical Test 2 ---------------------------------------------------

    % analytical test for the finally-operator

    % dynamic system, parameters, and options
    sys = linearSys([0 -1; 1 0],[0;0]);

    params.R0 = zonotope([0;-1]);
    params.tFinal = 2.5;

    options.timeStep = 0.5;
    options.zonotopeOrder = 10;
    options.taylorTerms = 10;

    % STL formula
    y = stl('y',2);
    goal = halfspace([0, -1],-0.4);
    eq = finally(in(y,goal),interval(1,2));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if ~modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end


    % Analytical Test 3 ---------------------------------------------------

    % analytical test for the globally-operator

    % dynamic system, parameters, and options
    sys = linearSys([0 -1; 1 0],[0;0]);

    params.R0 = zonotope([0;-1]);
    params.tFinal = 2;

    options.timeStep = 0.5;
    options.zonotopeOrder = 10;
    options.taylorTerms = 10;

    % STL formula
    y = stl('y',2);
    goal = halfspace([0, 1],0.5);
    eq = globally(in(y,goal),interval(1,2));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if ~modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end


    % Analytical Test 4 ---------------------------------------------------

    % analytical test for the globally-operator

    A = [0.9 0; 0 -0.9];
    B = 1;
    c = [0; 4.8];
    sys = linearSys(A,B,c);

    params.R0 = zonotope(interval([1;-1],[2;1]));
    params.tFinal = 3;

    options = struct();
    options.linAlg = 'adaptive';
    options.error = 0.1;

    % STL formula
    x = stl('x',2);
    eq = globally(x(2) < 5, interval(0,3));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end


    % Analytical Test 5 ---------------------------------------------------

    % analytical test for globally operator

    params.R0 = zonotope(interval([0.9;-0.1],[1.1;0.0]));
    params.tFinal = 1;

    options = struct();
    options.linAlg = 'adaptive';
    options.error = 0.1;

    A = diag([1 2]);
    B = [1;1];

    sys = linearSys(A,B);

    % STL formula
    x = stl('x',2);
    eq = globally(x(2) < 0.5,interval(0,1));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if ~modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end


    % Analytical Test 6 ---------------------------------------------------

    % analytical test for globally operator

    % hybrid system with two branches
    dyn = linearSys([0 0; 0 0], 0, [1;1]);

    reset.A = diag([1 1]);
    reset.c = [0;0];

    guard = conHyperplane([1 0],5);
    trans(1) = transition(guard,reset,2);

    guard = conHyperplane([0 1],5);
    trans(2) = transition(guard,reset,3);

    inv = polytope(diag([1 1]),[5;5]);
    loc(1) = location(inv,trans,dyn);

    dyn = linearSys(diag([0 0]), 0, [1;0]);
    inv = halfspace([-1 0], -5);
    loc(2) = location(inv,transition(),dyn);

    dyn = linearSys(diag([0 0]), 0, [0;1]);
    inv = halfspace([0 -1], -5);
    loc(3) = location(inv,transition(),dyn);

    sys = hybridAutomaton(loc);

    params.tFinal = 10;
    params.startLoc = 1;
    params.R0 = zonotope(interval([-0.1; -0.1], [0.1; 0.1]));

    options = struct();
    options.timeStep = 0.1;
    options.taylorTerms = 10;
    options.zonotopeOrder = 20;

    options.guardIntersect = 'polytope';
    options.enclose = {'box'};

    % STL formula
    x = stl('x',2);
    eq = globally(x(1) < 6, interval(0,10)) | ....
        globally(x(2) < 6, interval(0,10));

    % reachability analysis
    R = reach(sys,params,options);

    % model checking
    for i = 1:length(alg)
        if ~strcmp(alg{i}, 'signals')
            % 'rtl' and 'sampledTime' cannot handle branching right now
            continue
        end

        if ~modelChecking(R,eq,alg{i})
            throw(CORAerror('CORA:testFailed'));
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
