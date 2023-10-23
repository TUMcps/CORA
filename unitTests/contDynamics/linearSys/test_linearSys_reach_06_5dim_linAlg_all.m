function res = test_linearSys_reach_06_5dim_linAlg_all()
% test_linearSys_reach_06_5dim_linAlg_all - unit_test_function of linear reachability 
%    analysis with uncertain inputs (toy example),
%    all linear reach algorithms used (except krylov due to system size)
%    reduced version of testLong_linearSys_reach_06
%
%
% Syntax:
%    res = test_linearSys_reach_06_5dim_linAlg_all()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       26-June-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
dim_x = length(A);
fiveDimSys = linearSys('fiveDimSys',A,B);


% Parameters --------------------------------------------------------------

params.tFinal = 1;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = 0.5*zonotope([1; 0; 0; 0.5; -0.5],diag([0.2, 0.5, 0.2, 0.5, 0.5]));


% Reachability Settings ---------------------------------------------------

options.timeStep      = 0.05;
options.taylorTerms   = 4; % has to be 3 (calculated in reach_decomp)
options.zonotopeOrder = 20;

spec = specification(halfspace([0 -1 0 0 0],-2));


% Reachability Analysis ---------------------------------------------------

algs = {'standard','wrapping-free','fromStart','decomp'};
% note: krylov does not make sense for this small system dimension
Rs = cell(length(algs),1);
intHull = cell(length(algs),1);

for i=1:length(algs)
	options.linAlg = algs{i};
    if strcmp(options.linAlg,'decomp')
        % additional options when decomp called
        options.partition = [1, 2; 3, 4; 5, 5];
        Rs{i} = reach(fiveDimSys, params, options, spec);
        intHull{i} = interval(Rs{i}.timeInterval.set{end});
        options = rmfield(options,'partition');
    else % 'standard', 'wrapping-free', 'fromStart'
        Rs{i} = reach(fiveDimSys, params, options, spec);
        intHull{i} = interval(Rs{i}.timeInterval.set{end});
    end
end


% Simulation --------------------------------------------------------------

simOpt.points = 5;
simRes = simulateRandom(fiveDimSys, params, simOpt);


% Visualization -----------------------------------------------------------

plotting = false;
overlap  = false;
colors = {'b';'r';'g';'m'};
projDims = [1,2];

if plotting
    
    for i=1:length(algs)
        
        if ~overlap
            figure;
        else
            hold on
        end
        
        % plot initial set
        plot(params.R0,projDims,'b-','lineWidth',2);
        hold on

        % plot reach results
        if ~strcmp(algs{i},'decomp')
            plot(Rs{i},projDims,'FaceColor',[.9 .9 .9],'EdgeColor',colors{i});
        else
            % decomp
            % dimensions [1,2] are in (1)
            plot(Rs{i}(1),projDims,'FaceColor',[.9 .9 .9],'EdgeColor',colors{i});
        end
        
        % plot simulation results
        if ~overlap
            plot(simRes,projDims,'k');
            % label plot
            xlabel(['x_{',num2str(projDims(1)),'}']);
            ylabel(['x_{',num2str(projDims(2)),'}']);
        end
        
    end
    
    if overlap
        plot(simRes,projDims,'k');
        xlabel(['x_{',num2str(projDims(1)),'}']);
        ylabel(['x_{',num2str(projDims(2)),'}']);
    end
    
    hold off
end


% Numerical Comparison ----------------------------------------------------

% notes:
% - wrapping-free overapproximates standard if inhomogeneous solution
% - fromStart and decomp have to be exactly equal
% - standard and fromStart will be different because standard reduces
%    hom. and part. solution, whereas fromStart reduces only part. solution

% order should be: standard, wrapping-free, fromStart, decomp
% coloring: blue, red, green, magenta

plotHulls = false;
if plotHulls
    figure;
    plot(intHull{1},projDims,colors{1});
    hold on
    for i=2:length(algs)
        plot(intHull{i},projDims,colors{i});
    end
end

% check first one with all others
for i=2:length(algs)
    % check if slightly bloated versions enclose each other
    res_1 = (intHull{1} <= enlarge(intHull{i},1+1e-8));
    res_2 = (intHull{i} <= enlarge(intHull{1},1+1e-8)); % 0 for wrapping-free

    % final result
    if i == 2
        % wrapping-free contains standard, but not vice versa
        res_zono = res_1;
    else
        % standard - fromStart only differ due to reduction (within bloat)
        % fromStart - decomp should be equal
        res_zono = res_1 && res_2;
    end
    if ~res_zono
        res = false;
        break
    end
end


% ------------------------------ END OF CODE ------------------------------
