function res = example_nonlinear_reach_14_adaptiveHSCC2
% example_nonlinear_reach_14_adaptiveHSCC2 - example for nonlinear
%    reachability analysis using adaptive parameter tuning,
%    reproducing results from [1]
%
% Syntax:
%    res = example_nonlinear_reach_14_adaptiveHSCC2
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, M. Althoff. "Adaptive Parameter Tuning
%        for Reachability Analysis of Nonlinear Systems", HSCC 2021.

% Authors:       Mark Wetzlinger
% Written:       02-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% system and options from subfiles
prefix = 'aux_adaptive_';
handles = [ % dimensions
    "jetEngine()"; % 2
    "vanDerPol()"; % 2
    "brusselator()"; % 2
    "roessler()"; % 3
    "lorenz()"; % 3
    "springpendulum()"; % 4
    "lotkaVolterra()"; % 5
    "biologicalModel()"; % 7
    "genetic()"; % 9
    ];
sys_total = length(handles);

% output variable: tab3
tab3 = zeros(length(handles),12);
% columns (6 for lin + 6 for poly)
timeIdx = 1;
deltatminIdx = 2;
deltatmaxIdx = 3;
rhomaxIdx = 4;
lmaxIdx = 5;
gammaminIdx = 6;


% loop over all nonlinear benchmarks
for sys_no = 1:sys_total
    [sys, params] = eval(prefix + handles(sys_no));
    disp("System " + sys_no + ": " + sys.name);
    algs = {'lin-adaptive','poly-adaptive'};

    for a=1:length(algs)

        options.alg = algs{a};
        options.verbose = true;

        adapTime = tic;
        [R,~,opt] = reach(sys,params,options);
        endset = R.timePoint.set{end};
        tVec = query(R,'tVec');

        tComp = toc(adapTime);

        % simulation
        simOpt.points = 1000;                       % number of initial points
        simOpt.fracVert = 2^sys.dim/simOpt.points;  % fraction of vertices initial set

        simRes = simulateRandom(sys,params,simOpt);

        % computation of gamma_min
        endpoints = zeros(sys.dim,simOpt.points);
        for i=1:simOpt.points
            endpoints(:,i) = simRes(i).x{1}(end,:)';
        end
        simendset = interval.enclosePoints(endpoints);
        if contains(options.alg,'poly')
            methods = {'interval'}; %,'bnb','bnbAdv','globOpt'};
            edgelengths = [];
            for i=1:length(methods)
                try
                    edgelengths(:,end+1) = 2*rad(interval(endset,methods(i)));
                catch
                    edgelengths(:,end+1) = Inf(sys.dim,1);
                end
            end
            [gamma_o,minidx] = min(edgelengths,[],2);
        else
            gamma_o = 2*rad(interval(endset));
        end

        gamma_u = 2*rad(interval(simendset));
        gamma_min = min(gamma_u ./ gamma_o);


        % metrics for table -------------------------------------------
        tab3(sys_no,6*(a-1)+timeIdx) = tComp;
        tab3(sys_no,6*(a-1)+deltatminIdx) = min(tVec);
        tab3(sys_no,6*(a-1)+deltatmaxIdx) = max(tVec);
        tab3(sys_no,6*(a-1)+rhomaxIdx) = max(sum(opt.zonordersRtp,2));
        tab3(sys_no,6*(a-1)+lmaxIdx) = max(gamma_o);
        tab3(sys_no,6*(a-1)+gammaminIdx) = gamma_min;
        % -------------------------------------------------------------

    end
    
end

end


% Auxiliary functions -----------------------------------------------------

% Investigated Systems ----------------------------------------------------

function [sys, params] = aux_adaptive_vanDerPol()

params.tFinal = 6.74;
params.R0 = zonotope([[1.4;2.4],diag([0.14,0.05])]);
sys = nonlinearSys(@vanderPolEq,2,1);

end

function [sys, params] = aux_adaptive_brusselator()

params.tFinal = 5; 
params.R0 = zonotope(interval([0.9;0],[1;0.1]));
sys = nonlinearSys(@brusselator,2,1);

end

function [sys, params] = aux_adaptive_jetEngine()

params.tFinal = 8;
params.R0 = zonotope([[1;1],0.1*diag(ones(2,1))]);
sys = nonlinearSys(@jetEngine,2,1);

end

function [sys, params] = aux_adaptive_lorenz()

params.tFinal = 2;
params.R0 = zonotope([[15;15;35],0.1*diag(ones(3,1))]);
sys = nonlinearSys(@lorenz,3,1);

end

function [sys, params] = aux_adaptive_biologicalModel()

params.tFinal = 2;
params.R0 = zonotope([ones(7,1),0.01*diag(ones(7,1))]);
sys = nonlinearSys(@biologicalModel,7,1);

end

function [sys, params] = aux_adaptive_roessler()

params.tFinal = 6;
params.R0 = zonotope([[0;-8.4;0],0.2*diag(ones(3,1))]);
sys = nonlinearSys(@roessler,3,1);

end

function [sys, params] = aux_adaptive_lotkaVolterra()

params.tFinal = 5;
params.R0 = zonotope([0.95*ones(5,1),0.05*diag(ones(5,1))]);
sys = nonlinearSys(@lotkaVolterraCont,5,1);

end

function [sys, params] = aux_adaptive_genetic()

params.tFinal = 0.1;
centerR0 = [1;1.3;0.1;0.1;0.1;1.3;2.5;0.6;1.3];
W0 = 0.04; % 0.01 | 0.02 | 0.04
params.R0 = zonotope(centerR0,diag(W0*ones(9,1)));
sys = nonlinearSys(@genetic,9,1);

end

function [sys, params] = aux_adaptive_springpendulum()

params.tFinal = 1;
params.R0 = zonotope([1.2;0.5;0.05;0.05],diag([0.1,0.1,0.05,0.05]));
sys = nonlinearSys(@springpendulum,4,1);

end

% ------------------------------ END OF CODE ------------------------------
