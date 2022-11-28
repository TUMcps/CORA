function completed = run_arch22_ainncs_category()
% run_arch22_ainncs_category - runs all benchmark of the
%                              ARCH'22 AINNCS Category
%
% Syntax:
%    completed = run_arch22_ainncs_category()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP22 Category Report:
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants."
%       EPiC Series in Computing 80 (2022): 90-119.
%       TODO: correct reference

% Author:       Tobias Ladner
% Written:      20-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% STARTUP -----------------------------------------------------------------
rng(1)
warning off

basepath = "./results";
mkdir(basepath)
plotpath = sprintf("%s/plots", basepath);
mkdir(plotpath)

% for smooth images
set(0, 'defaultFigureRenderer', 'painters')

resultstxt = sprintf("%s/results.txt", basepath);
delete(resultstxt)
diary(resultstxt)
disp("--------------------------------------------------------")
disp("ARCH'22 AINNCS Category")
disp("Tool: CORA - A Tool for Continuous Reachability Analysis")
fprintf("Date: %s\n", datestr(datetime()))
disp("--------------------------------------------------------")
disp(" ")

% RUN BENCHMARKS ----------------------------------------------------------
benchmarks = {; ...
    @example_neuralNet_reach_01_unicycle, "01_unicycle"; ...
    @example_neuralNet_reach_02_ACC, "02_ACC"; ...
    @example_neuralNet_reach_03_TORA, "03_TORA"; ...
    @example_neuralNet_reach_03_TORA_heterogeneousTanh, "03_TORA_heterogeneousTanh"; ...
    @example_neuralNet_reach_03_TORA_heterogeneousSigmoid, "03_TORA_heterogeneousSigmoid"; ...
    @example_neuralNet_reach_04_singlePendulum, "04_singlePendulum"; ...
    @example_neuralNet_reach_05_doublePendulum_lessRobust, "05_doublePendulum_lessRobust"; ...
    @example_neuralNet_reach_05_doublePendulum_moreRobust, "05_doublePendulum_moreRobust"; ...
    @example_neuralNet_reach_06_airplane, "06_airplane"; ...
    @example_neuralNet_reach_07_VCAS_middle, "07_VCAS_middle"; ...
    @example_neuralNet_reach_07_VCAS_worst, "07_VCAS_worst"; ...
    @example_neuralNet_reach_08_attitudeControl, "08_attitudeControl"; ...
    };

n = size(benchmarks, 1);
fprintf("Running %d Benchmarks.. \n", n);
disp(" ")

for i = 1:n
    disp("--------------------------------------------------------")
    benchmark = benchmarks{i, 1};
    name = benchmarks{i, 2};

    try
        % run benchmark
        benchmark();

        disp(" ")
        fprintf("'%s' was run successfully!\n", name)
        fprintf("Saving plots to '%s'..\n", plotpath)

        % save plots
        h = findobj('type', 'figure');
        m = length(h);

        for j = 1:m
            saveas(gcf, sprintf("./%s/cora_%s_%d.%s", plotpath, name, j, 'png'));
            saveas(gcf, sprintf("./%s/cora_%s_%d.%s", plotpath, name, j, 'eps'), 'epsc');
            close(gcf)
        end

    catch ME
        % error handling
        disp(" ")
        fprintf("An ERROR occured during execution of '%s':\n", name);
        disp(ME.getReport())
        disp("Continuing with next benchmark..")
    end

    disp(" ")
end

% -------------------------------------------------------------------------
disp("--------------------------------------------------------")
disp(" ")
completed = 1;
disp("Completed!")
diary off;
end

%------------- END OF CODE --------------