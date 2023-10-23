function res = test_stl_modelCheckTrace()
% test_stl_modelCheckTrace - unit test function of modelCheckTrace
%
% Syntax:
%    res = test_stl_modelCheckTrace
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

    % Analytical Test 1 ---------------------------------------------------

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
    if modelCheckTrace(eq,x,t) ~= true
        throw(CORAerror('CORA:testFailed'));
    end

    % slightly change STL formula
    y = stl('y',2);
    point = sqrt(2)/2 + 0.01;
    eq = until(y(2) < -point,y(1) > point,interval(0,1));

    % model checking
    if modelCheckTrace(eq,x,t) ~= false
        throw(CORAerror('CORA:testFailed'));
    end


    % analytical tests for the globally-operator

    % trajectory
    x = linspace(0,1,100)';
    t = linspace(0,1,100)';

    % STL formula
    y = stl('y',2); 
    eq = globally(y(1) > 0.1,interval(0,1));

    % model checking
    if modelCheckTrace(eq,x,t) ~= false
        throw(CORAerror('CORA:testFailed'));
    end

    % different STL formula
    eq = globally(y(1) >= 0.2 & y(1) <= 0.3, interval(0.2,0.3));

    % model checking
    if modelCheckTrace(eq,x,t) ~= true
        throw(CORAerror('CORA:testFailed'));
    end   


    % analytical tests for the finally-operator

    % trajectory
    x = linspace(0,1,100)';
    t = linspace(0,1,100)';

    % STL formula
    y = stl('y',2); 
    eq = finally(y(1) > 1.1,interval(0,1));

    % model checking
    if modelCheckTrace(eq,x,t) ~= false
        throw(CORAerror('CORA:testFailed'));
    end

    % different STL formula
    eq = finally(y(1) >= 0.25, interval(0.2,0.3));

    % model checking
    if modelCheckTrace(eq,x,t) ~= true
        throw(CORAerror('CORA:testFailed'));
    end   


    % Random Tests --------------------------------------------------------

    for i = 1:10

        % generate random linear system
        sys = linearSys(-1+2*rand(2,2),-1+2*rand(2,1));
    
        % generate random trace
        simOpts.x0 = -2+4*rand(2,1);
        simOpts.u = -2 + 4*rand(1,10);
        simOpts.tFinal = 1;
    
        [t_,x_] = simulate(sys,simOpts);
    
        t = 0:0.01:1;
        [~,ind] = unique(t_);
        x = [interp1(t_(ind),x_(ind,1),t,'linear','extrap'); ...
             interp1(t_(ind),x_(ind,2),t,'linear','extrap')]';
        
        % generate random STL formula
        dom = interval(min(x,[],1)',max(x,[],1)');
        eq = stl.generateRandom('Dimension',2,'FinalTime',1, ...
                                            'TimeStep',0.1,'Domain',dom);
    
        % model check original formula
        resOrig = modelCheckTrace(eq,x,t);

        % model check inverted formula
        resInv = modelCheckTrace(~eq,x,t);

        % check if results are consistent
        if resOrig == resInv
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
