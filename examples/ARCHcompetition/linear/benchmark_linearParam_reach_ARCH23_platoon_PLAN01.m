function text = benchmark_linearParam_reach_ARCH23_platoon_PLAN01()
% benchmark_linearParam_reach_ARCH23_platoon_PLAN01 - example of linear 
%    reachability analysis from the ARCH23 friendly competition
%    (platoon example); the linear dynamics switch nondeterministically
%
% Syntax:
%    benchmark_linearParam_reach_ARCH23_platoon_PLAN01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       09-February-2017
% Last update:   13-March-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameter ---------------------------------------------------------------

% initial set
R0 = zonotope(zeros(9,1));
params.R0 = R0;

% uncertain inputs
params.U = zonotope(interval(-9,1));

% final time
params.tFinal = 50;


% Reachability Settings ---------------------------------------------------

options.intermediateTerms = 3;
options.timeStep = 0.01;
options.zonotopeOrder = 400;
options.taylorTerms = 3;
options.compTimePoint = false;


% System Dynamics ---------------------------------------------------------

% controlled system
A_c = [...
        0    1.0000         0         0         0         0         0         0         0;...
        0         0   -1.0000         0         0         0         0         0         0;...
        1.6050    4.8680   -3.5754   -0.8198    0.4270   -0.0450   -0.1942    0.3626   -0.0946;...
        0         0         0         0    1.0000         0         0         0         0;...
        0         0    1.0000         0         0   -1.0000         0         0         0;...
        0.8718    3.8140   -0.0754    1.1936    3.6258   -3.2396   -0.5950    0.1294   -0.0796;...
        0         0         0         0         0         0         0    1.0000         0;...
        0         0         0         0         0    1.0000         0         0   -1.0000;...
        0.7132    3.5730   -0.0964    0.8472    3.2568   -0.0876    1.2726    3.0720   -3.1356 ]; 

% uncontrolled system
A_n = [...
        0    1.0000         0         0         0         0         0         0         0;...
        0         0   -1.0000         0         0         0         0         0         0;...
        1.6050    4.8680   -3.5754         0         0         0         0         0         0;...
        0         0         0         0    1.0000         0         0         0         0;...
        0         0    1.0000         0         0   -1.0000         0         0         0;...
        0         0         0    1.1936    3.6258   -3.2396         0         0         0;...
        0         0         0         0         0         0         0    1.0000         0;...
        0         0         0         0         0    1.0000         0         0   -1.0000;...
        0.7132    3.5730   -0.0964    0.8472    3.2568   -0.0876    1.2726    3.0720   -3.1356 ];   
    
    
B = [0 ; 1; 0; 0; 0; 0; 0; 0; 0 ];
        
% build zonotpe matrix
A_mid = 0.5*(A_c + A_n);
A_gen{1} = 0.5*(A_c - A_n);
matZ_A = matZonotope(A_mid, A_gen);

% instantiate linear dynamics with constant parameters
linSys  = linParamSys(matZ_A, B,'varParam');


% Reachability Analysis ---------------------------------------------------

% reachable set computations
tic
R = reach(linSys, params, options);

% invariant computation
Rlast = reduce(R.timeInterval.set{end},'girard', 400);
options.R0 = enlarge(Rlast,1.01);
options.zonotopeOrder = 800;
options.U = params.U;

[RcontInv,tadd] = aux_reachInv(linSys, options, options.R0);

tComp = toc;


% Verification ------------------------------------------------------------

tic
violation50 = false;
for i = 1:length(R.timeInterval.set)
    x_proj = interval(project(R.timeInterval.set{i},[1,4,7]));
    if any(infimum(x_proj) < -50)
        violation50 = true;
    end
end
for i=1:length(RcontInv)
    x_proj = interval(project(RcontInv{i},[1,4,7]));
    if any(infimum(x_proj) < -50)
        violation50 = true;
    end
end

tVer = toc;

disp(['specifications verified: ',num2str(~violation50)]);
disp(['computation time: ',num2str(tVer + tComp)]);


% Simulation --------------------------------------------------------------

% update parameter
params.R0 = R0;
params.tFinal = params.tFinal + tadd;

% simulation options
simOpt.points = 10;

% random simulation
simRes = simulateRandom(linSys, params, simOpt);


% Visualization -----------------------------------------------------------

figure; hold on; box on;

useCORAcolors('CORA:contDynamics');

% plot reachable set over time (part 1)
plotOverTime(R,1,'Unify',true);

RcontInvOverTime = cartProd(interval(params.tFinal-options.timeStep),...
    interval(project(RcontInv{i},1)));
plot(RcontInvOverTime,1,'FaceColor',CORAcolor("CORA:reachSet"));

% plot simulation results
plotOverTime(simRes,1);
% formatting
xlabel('$t$','interpreter','latex');
ylabel('$x_1$','interpreter','latex');
xlim([0,50]);

% return value
text = ['Platoon,PLAN01-UNB50,',num2str(~violation50),',',num2str(tVer + tComp)];

end


% Auxiliary functions -----------------------------------------------------

function [Rcont,t] = aux_reachInv(obj,options,inv)
% compute the reachable set until it is fully contained inside the
% invariant set

    % adapt options
    options.uTrans = center(options.U);
    options.U = options.U - options.uTrans;
    options.originContained = 1;
    options.reductionTechnique = 'girard';
    options.tStart = 0;

    % obtain factors for initial state and input solution
    for i = 1:(options.taylorTerms+1)
        r = options.timeStep;
        options.factor(i)= r^(i)/factorial(i);    
    end

    % initialize reachable set computations
    [Rnext, options] = initReach(obj, options.R0, options);

    % loop until not in invariant set anymore
    t = options.tStart;
    iSet = 1;
    notInInv = 1;

    while notInInv

        % save reachable set in cell structure
        Rcont{iSet} = Rnext.ti; 
        Rcont_tp{iSet} = Rnext.tp; 

        % increment time and set counter
        t = t + options.timeStep;
        iSet = iSet + 1; 
        options.t = t;

        % compute next reachable set
        [Rnext,options]=post(obj,Rnext,options);

        % check if reachable set is contained in invariant
        if aux_inViaProj(Rnext.ti, inv)
            notInInv=0;
        end
    end
end

function res = aux_inViaProj(Z1,Z2)
% check if zonotoep Z1 is contained in Z2

    % dimension
    dimension = length(center(Z1));

    % combinations of projections
    comb = combinator(dimension,2,'c');

    % init result
    res = true;

    % loop through projections
    for proj = 1:length(comb(:,1))
        
        % perform projections
        Z1proj = project(Z1,comb(proj,:));
        Z2proj = project(Z2,comb(proj,:));

        % check enclosure
        if ~contains(Z1proj, Z2proj)
            res = false;
            break
        end
    end  
end

% ------------------------------ END OF CODE ------------------------------
