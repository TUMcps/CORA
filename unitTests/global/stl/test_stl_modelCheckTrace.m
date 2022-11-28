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
%    res - boolean 
%
% Example: 
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
        error('Analytical test 1 failed!')
    end

    % slightly change STL formula
    y = stl('y',2);
    point = sqrt(2)/2 + 0.01;
    eq = until(y(2) < -point,y(1) > point,interval(0,1));

    % model checking
    if modelCheckTrace(eq,x,t) ~= false
        error('Analytical test 1 failed!')
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

            % save variables so that failure can be reproduced
            path = pathFailedTests(mfilename());
            save(path,'eq','x','t');

            error('Random test failed!');
        end
    end
end

%------------- END OF CODE --------------