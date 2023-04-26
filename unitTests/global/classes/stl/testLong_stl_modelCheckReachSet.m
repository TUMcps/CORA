function res = testLong_stl_modelCheckReachSet()
% testLong_stl_modelCheckReachSet - unit test function of
%    modelCheckReachSet
%
% Syntax:  
%    res = testLong_stl_modelCheckReachSet
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

% Analytical Test ---------------------------------------------------------

% analytical test for the until-operator

% dynamic system, parameters, and options
sys = linearSys([0 -1; 1 0],[0;0]);

params.R0 = zonotope([0;-1]);
params.tFinal = 2;

options.timeStep = 1;
options.zonotopeOrder = 10;
options.taylorTerms = 10;

% STL formula
y = stl('y',2);
point = sqrt(2)/2;
eq = until(y(2) < -point,y(1) > point,interval(0,2));

% loop over different time step sizes
for i = 1:4

    % reachability analysis
    options.timeStep = options.timeStep/2;

    R = reach(sys,params,options);

    % model checking (should be false for arbitrary time steps)
    if modelCheckReachSet(eq,R) ~= false
        throw(CORAerror('CORA:testFailed'));
    end
end

% modified STL formula
y = stl('y',2);
point = 0.7;
eq = until(y(2) < -point,y(1) > point,interval(0,2));

% refine time step until equation satisfied
options.timeStep = 1; resTmp = false;

for i = 1:5

    % reachability analysis
    options.timeStep = options.timeStep/2;

    R = reach(sys,params,options);

    % model checking (should be false for arbitrary time steps)
    if modelCheckReachSet(eq,R) == true
        resTmp = true; break;
    end
end

if ~resTmp
    throw(CORAerror('CORA:testFailed'));
end


% Random Tests ------------------------------------------------------------

cnt = 0;

while cnt < 5

    % generate random linear system
    sys = linearSys(-1+2*rand(2,2),-1+2*rand(2,1));

    % generate random trace
    simOpts.x0 = -2+4*rand(2,1);
    simOpts.u = -2 + 4*rand(1,5);
    simOpts.tFinal = 1;

    [t_,x_] = simulate(sys,simOpts);

    t = 0:0.01:1;
    [~,ind] = unique(t_);
    x = [interp1(t_(ind),x_(ind,1),t,'linear','extrap'); ...
         interp1(t_(ind),x_(ind,2),t,'linear','extrap')]';
    
    % generate random STL formula
    dom = interval(min(x,[],1)',max(x,[],1)');
    nrOps = randi(5);
    eq = stl.generateRandom('Dimension',2,'FinalTime',1,'Domain',dom, ...
                            'NrOperators',nrOps,'TimeStep',0.2,'NestedOps',1);

    % model check the trace
    resTrace = modelCheckTrace(eq,x,t);

    % if the trace violates the formula, the reachable set also has to
    if ~resTrace

        cnt = cnt + 1;

        % compute reachable set
        params.R0 = zonotope(simOpts.x0);
        params.u = simOpts.u;
        params.tFinal = simOpts.tFinal;

        options.timeStep = 0.2;
        options.zonotopeOrder = 10;
        options.taylorTerms = 4;

        R = reach(sys,params,options);

        % model check reachable set
        resReach = modelCheckReachSet(eq,R);

        % check if results are consistent
        if resReach
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

%------------- END OF CODE --------------