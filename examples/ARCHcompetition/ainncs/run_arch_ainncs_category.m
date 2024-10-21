function completed = run_arch_ainncs_category()
% run_arch_ainncs_category - runs all benchmark of the ARCH AINNCS Category
%
% Syntax:
%    completed = run_arch_ainncs_category()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Reference:
%   [1] Lopez, Diego M., et al. "ARCH-COMP22 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants." Proceedings of 
%       9th International Workshop on Applied. Vol. 90. 2022.

% Authors:       Tobias Ladner
% Written:       20-May-2022
% Last update:   30-March-2023 (TL, ARCH'23 revisions)
%                02-May-2024 (TL, ARCH'24 revisions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% STARTUP -----------------------------------------------------------------
rng(1)
warning off

basepath = "./results";
mkdir(basepath)
plotpath = sprintf("%s/plots", basepath);
mkdir(plotpath)

% for smooth images
set(0, 'defaultFigureRenderer', 'painters')

% diary
resultstxt = sprintf("%s/results.txt", basepath);
delete(resultstxt)
diary(resultstxt)

disp("--------------------------------------------------------")
disp("ARCH'24 AINNCS Category")
disp("Tool: CORA - A Tool for Continuous Reachability Analysis")
fprintf("Date: %s\n", datestr(datetime()))
disp("--------------------------------------------------------")
disp(" ")

% RUN BENCHMARKS ----------------------------------------------------------
% each benchmark returns [completed, res, tTotal]
% and might visualize several plots

benchmarks = {; ...
    % {@func, 'benchmark', 'instance'}
    @example_neuralNet_reach_01_unicycle, 'Unicycle', 'reach'; ...
    @example_neuralNet_reach_02_ACC, 'ACC', 'safe-distance'; ...
    @example_neuralNet_reach_03_TORA, 'TORA', 'remain'; ...
    @example_neuralNet_reach_03_TORA_heterogeneousTanh, 'TORA', 'reach-tanh'; ...
    @example_neuralNet_reach_03_TORA_heterogeneousSigmoid, 'TORA', 'reach-sigmoid'; ...
    @example_neuralNet_reach_04_singlePendulum, 'Single_Pendulum', 'reach'; ...
    @example_neuralNet_reach_05_doublePendulum_lessRobust, 'Double_Pendulum', 'less-robust'; ...
    @example_neuralNet_reach_05_doublePendulum_moreRobust, 'Double_Pendulum', 'more-robust'; ...
    @example_neuralNet_reach_06_airplane, 'Airplane', 'continuous'; ...
    @() example_neuralNet_reach_07_VCAS('middle', -19.5), 'VCAS', 'middle-19.5'; ...
    @() example_neuralNet_reach_07_VCAS('middle', -22.5), 'VCAS', 'middle-22.5'; ...
    @() example_neuralNet_reach_07_VCAS('middle', -25.5), 'VCAS', 'middle-25.5'; ...
    @() example_neuralNet_reach_07_VCAS('middle', -28.5), 'VCAS', 'middle-28.5'; ...
    @() example_neuralNet_reach_07_VCAS('worst', -19.5), 'VCAS', 'worst-19.5'; ...
    @() example_neuralNet_reach_07_VCAS('worst', -22.5), 'VCAS', 'worst-22.5'; ...
    @() example_neuralNet_reach_07_VCAS('worst', -25.5), 'VCAS', 'worst-25.5'; ...
    @() example_neuralNet_reach_07_VCAS('worst', -28.5), 'VCAS', 'worst-28.5'; ...
    @example_neuralNet_reach_08_attitudeControl, 'Attitude_Control', 'avoid'; ...
    @benchmark_neuralNet_reach_09_QUAD, 'QUAD', 'reach'; ...
    @example_neuralNet_reach_10_spacecraftDocking, 'Docking', 'constraint'; ...
    @() example_neuralNet_reach_11_NAV('point'), 'NAV', 'standard'; ...
    @() example_neuralNet_reach_11_NAV('set'), 'NAV', 'robust'; ...
    @example_neuralNet_reach_12_cartPole, 'CartPole*', 'reach'; ...
    };

n = size(benchmarks, 1);
fprintf("Running %d Benchmarks.. \n", n);
disp(" ")

% results.csv
resultsCSV = cell(size(benchmarks, 1)+1, 4);
resultsCSV(1, :) = {'benchmark','instance','result','time'};

for i = 1:n
    disp("--------------------------------------------------------")
    benchmark = benchmarks{i, 1};
    name = benchmarks{i, 2};
    instance = benchmarks{i, 3};

    benchID = [num2str(i,'%02d') '_' name '_' instance];
    if isempty(instance)
        benchID = benchID(1:end-1); % remove trailing '_'
    end

    try
        % run benchmark
        [~, res, tTotal] = benchmark();
        resultsCSV(i+1, :) = {name, instance, res, num2str(tTotal)};

        disp(" ")
        fprintf("'%s' was run successfully!\n", benchID)
        fprintf("Saving plots to '%s'..\n", plotpath)

        % save plots
        h = findobj('type', 'figure');
        m = length(h);

        for j = 1:m
            plotFile = sprintf("./%s/cora_%s_%02d", plotpath, strrep(benchID,'*',''), j);
            savefig(gcf, sprintf("%s.%s", plotFile, 'fig'))
            saveas(gcf, sprintf("%s.%s", plotFile, 'png'))
            % saveas(gcf, sprintf("%s.%s", plotFile, 'eps'), 'epsc')
            close(gcf)
        end

    catch ME
        % error handling
        disp(" ")
        fprintf("An ERROR occured during execution of '%s':\n", benchID);
        disp(ME.getReport())
        disp("Continuing with next benchmark..")
    end

    disp(" ")
end

% write results.csv
writecell(resultsCSV, sprintf("%s/results.csv", basepath))

% -------------------------------------------------------------------------

disp("--------------------------------------------------------")
disp(" ")
completed = 1;
disp("Completed!")
diary off;
end

% ------------------------------ END OF CODE ------------------------------
