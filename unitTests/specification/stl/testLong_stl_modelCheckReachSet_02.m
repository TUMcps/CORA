function res = testLong_stl_modelCheckReachSet_02()
% testLong_stl_modelCheckReachSet_02 - unit test function of
%    modelCheckReachSet
%
% Syntax:
%    res = testLong_stl_modelCheckReachSet_02
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

% assume true
res = true;

% algorithms for reachSet model checking
alg = {'rtl','sampledTime','signals'};

% init counter
cnt = 0;

% check five cases
while cnt < 5

    % generate random 2D-linear system
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
        for j = 1:length(alg)

            resReach = modelChecking(R,eq,alg{j});
    
            % check if results are consistent
            if resReach
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
