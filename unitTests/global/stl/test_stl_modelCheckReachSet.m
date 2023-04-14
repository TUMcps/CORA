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

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    res = true;

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
    if modelCheckReachSet(eq,R) ~= false
        throw(CORAerror('CORA:testFailed'));
    end


    % Analytical Test 2 ---------------------------------------------------

    % analytical test for the finally-operator

    % dynamic system, parameters, and options
    sys = linearSys([0 -1; 1 0],[0;0]);
 
    params.R0 = zonotope([0;-1]);
    params.tFinal = 2;

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
    if modelCheckReachSet(eq,R) ~= true
        throw(CORAerror('CORA:testFailed'));
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
    if modelCheckReachSet(eq,R) ~= true
        throw(CORAerror('CORA:testFailed'));
    end

end

%------------- END OF CODE --------------