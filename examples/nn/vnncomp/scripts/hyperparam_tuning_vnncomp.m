function hyperparam_tuning_vnncomp(benchmarks,paramGrid,instanceIds,resumePath,datapath,timeoutMultiplier)
% hyperparam_tuning_vnncomp - grid search over hyperparameter values.
%
% Syntax:
%    hyperparam_tuning_vnncomp(benchmarks,paramGrid)
%    hyperparam_tuning_vnncomp(benchmarks,paramGrid,instanceIds)
%    hyperparam_tuning_vnncomp(benchmarks,paramGrid,instanceIds,resumePath)
%    hyperparam_tuning_vnncomp(benchmarks,paramGrid,instanceIds,resumePath,datapath)
%    hyperparam_tuning_vnncomp(benchmarks,paramGrid,instanceIds,resumePath,datapath,timeoutMultiplier)
%
% Inputs:
%    benchmarks - cell array of benchmark names
%    paramGrid - struct where each field is a param name and each value
%                is an array of values to try. Field names may use '__' to
%                address nested fields, e.g. 'train__mini_batch_size'.
%                Example: struct('num_splits',[2 4 8],'train__mini_batch_size',[16 64])
%                Each combination in the grid is the complete options.nn
%                config — no base defaults are merged in.
%    instanceIds - (optional) array of instance indices to run (default: all)
%    resumePath - (optional) path to a previous tuning results directory
%                 to resume from. Completed configs will be skipped.
%    datapath - (optional) path to the data directory (default: vnncomp/data)
%    timeoutMultiplier - (optional) fraction of the official per-instance
%                        timeout to use, e.g. 0.5 uses half the time
%                        (default: 1.0)
%
% Outputs:
%    -
%

% Authors:       Benedikt Kellner
% Written:       05-March-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Default: run all instances.
if nargin < 3 || isempty(instanceIds)
    instanceIds = [];
end
if nargin < 4
    resumePath = '';
end

if nargin < 6 || isempty(timeoutMultiplier)
    timeoutMultiplier = 1.0;
end

% Resolve paths relative to the vnncomp directory (parent of scripts/).
vnncompDir = fileparts(fileparts(mfilename('fullpath')));
if nargin < 5 || isempty(datapath)
    datapath = fullfile(vnncompDir,'data');
end

% Resume from previous run or create new results directory.
if ~isempty(resumePath)
    resultspath = resumePath;
else
    resultspath = fullfile(vnncompDir,'results',datestr(datetime,'yymmdd-hhMMss'),'tuning');
    mkdir(resultspath);
end

% Start diary to persist all command-window output to disk (append on resume).
logFile = fullfile(resultspath,'tuning.log');
diary(logFile);
cleanupDiary = onCleanup(@() diary('off'));

if ~isempty(resumePath)
    fprintf('Resuming from: %s\n',resultspath);
else
    fprintf('Results will be saved to: %s\n',resultspath);
end

% Save run metadata.
metadata = struct();
metadata.matlabVersion = version;
metadata.system = computer;
metadata.hasGPU = logical(exist('gpuDeviceCount','file')) && gpuDeviceCount > 0;
if metadata.hasGPU
    gpuInfo = gpuDevice;
    metadata.gpuName = gpuInfo.Name;
end
[~,gitHash] = system('git rev-parse --short HEAD');
metadata.gitHash = strtrim(gitHash);
metadata.paramGrid = paramGrid;
metadata.benchmarks = benchmarks;
metadata.instanceIds = instanceIds;
metadata.timestamp = datetime('now');
save(fullfile(resultspath,'metadata.mat'),'metadata');

% Extract parameter names and value arrays.
paramNames = fieldnames(paramGrid);
numParams = length(paramNames);
paramValues = cell(1,numParams);
for k = 1:numParams
    paramValues{k} = paramGrid.(paramNames{k});
end

% Compute Cartesian product of all parameter values.
% Convert each parameter's values to a cell array for uniform handling.
for k = 1:numParams
    if ~iscell(paramValues{k})
        paramValues{k} = num2cell(paramValues{k});
    end
end
% Build index grids over each parameter's values.
nVals = cellfun(@numel, paramValues);
idxArgs = arrayfun(@(n) 1:n, nVals, 'UniformOutput', false);
if numParams == 1
    idxGrids = {(1:nVals(1))'};
else
    [idxGrids{1:numParams}] = ndgrid(idxArgs{:});
end
numCombinations = numel(idxGrids{1});

% Store summary data.
configNames  = cell(numCombinations,1);
configLabels = cell(numCombinations,1);
numVerifAll = zeros(numCombinations,1);
numFalsAll = zeros(numCombinations,1);
numUnknownAll = zeros(numCombinations,1);
totalTimeAll = zeros(numCombinations,1);
% Store per-instance results for each config (cell of tables).
allResults = cell(numCombinations,1);

for c = 1:numCombinations
    % Use a short sequential name — full config is saved in tuning_results.mat.
    configName = sprintf('config_%0*d', ceil(log10(numCombinations+1)), c);
    configNames{c} = configName;
    subResultsPath = sprintf('%s/%s',resultspath,configName);

    % Build a human-readable label (for logging/plots, not used as path).
    labelParts = {};
    for k = 1:numParams
        val = paramValues{k}{idxGrids{k}(c)};
        if isnumeric(val) && ~isscalar(val)
            valStr = mat2str(val);
        else
            valStr = char(string(val));
        end
        labelParts{end+1} = sprintf('%s=%s',paramNames{k},valStr);
    end
    configLabels{c} = strjoin(labelParts,', ');

    fprintf('==================================================================\n');
    fprintf('CONFIG %d/%d: %s\n  %s\n',c,numCombinations,configName,configLabels{c});
    fprintf('==================================================================\n');

    % Skip if this config was already completed (resume support).
    if isfolder(subResultsPath) && isfile(fullfile(resultspath,'tuning_results.mat'))
        prev = load(fullfile(resultspath,'tuning_results.mat'),'tuningResults');
        if ismember(configName,prev.tuningResults.configNames)
            idx = find(strcmp(prev.tuningResults.configNames,configName),1);
            numVerifAll(c) = prev.tuningResults.numVerifAll(idx);
            numFalsAll(c) = prev.tuningResults.numFalsAll(idx);
            numUnknownAll(c) = prev.tuningResults.numUnknownAll(idx);
            totalTimeAll(c) = prev.tuningResults.totalTimeAll(idx);
            allResults{c} = prev.tuningResults.allResults{idx};
            fprintf('--- Skipping (already completed)\n');
            continue;
        end
    end

    % Build options purely from the paramGrid — no base defaults merged in.
    options = struct('nn',struct());
    for k = 1:numParams
        parts = strsplit(paramNames{k},'__');
        options.nn = setfield(options.nn,parts{:},paramValues{k}{idxGrids{k}(c)});
    end

    % Run benchmarks with fully resolved options — no merging in pipeline.
    run_benchmarks(benchmarks,datapath,subResultsPath,options,instanceIds,timeoutMultiplier);

    % Read results from generated CSV files.
    numVerif = 0; numFals = 0; numUnknown = 0; totalTime = 0;
    configResults = table();
    benchdirs = dir(subResultsPath);
    benchdirs = benchdirs([benchdirs.isdir] & ...
        ~ismember({benchdirs.name},{'.','..'}));
    for b = 1:length(benchdirs)
        csvFile = fullfile(benchdirs(b).folder,benchdirs(b).name,'results.csv');
        if isfile(csvFile)
            T = readtable(csvFile,'Delimiter',',');
            numVerif = numVerif + sum(strcmp(T.results,'unsat'));
            numFals = numFals + sum(strcmp(T.results,'sat'));
            numUnknown = numUnknown + sum(strcmp(T.results,'unknown'));
            totalTime = totalTime + sum(T.totalTimes);
            configResults = [configResults; T];
        end
    end
    numVerifAll(c) = numVerif;
    numFalsAll(c) = numFals;
    numUnknownAll(c) = numUnknown;
    totalTimeAll(c) = totalTime;
    allResults{c} = configResults;

    % Save intermediate results after each config (crash-safe).
    tuningResults.configNames = configNames(1:c);
    tuningResults.numVerifAll = numVerifAll(1:c);
    tuningResults.numFalsAll = numFalsAll(1:c);
    tuningResults.numUnknownAll = numUnknownAll(1:c);
    tuningResults.totalTimeAll = totalTimeAll(1:c);
    tuningResults.allResults = allResults(1:c);
    tuningResults.configLabels = configLabels(1:c);
    tuningResults.paramGrid = paramGrid;
    save(fullfile(resultspath,'tuning_results.mat'),'tuningResults');
end

% Print summary table.
fprintf('\n');
summaryTable = CORAtableParameters('Hyperparameter Tuning Summary');
summaryTable.printHeader();
for c = 1:numCombinations
    total = numVerifAll(c) + numFalsAll(c) + numUnknownAll(c);
    solved = numVerifAll(c) + numFalsAll(c);
    summaryTable.printContentRow(sprintf('%s  (%s)',configNames{c},configLabels{c}), ...
        sprintf('V:%d F:%d U:%d Solved:%d/%d [%.1f%%] Time:%.1fs', ...
        numVerifAll(c),numFalsAll(c),numUnknownAll(c), ...
        solved,total,solved/max(total,1)*100,totalTimeAll(c)));
    if c < numCombinations
        summaryTable.printMidBoundaryRow();
    end
end
summaryTable.printFooter();
fprintf('Results saved to: %s\n',resultspath);

% Plotting is best-effort: never abort a long sweep on a render failure.
try
plotspath = fullfile(resultspath,'plots');
mkdir(plotspath);

fontSize = 11;

% Figure 1: Stacked bar chart of verified/falsified/unknown per config.
f1 = figure('Name','Verification Results');
barData = [numVerifAll, numFalsAll, numUnknownAll];
b = bar(barData,'stacked');
b(1).FaceColor = [0.2 0.7 0.3]; % verified = green
b(2).FaceColor = [0.9 0.3 0.2]; % falsified = red
b(3).FaceColor = [0.6 0.6 0.6]; % unknown = gray
set(gca,'XTickLabel',configLabels,'XTickLabelRotation',45,'FontSize',fontSize);
ylabel('Number of Instances');
legend({'Verified','Falsified','Unknown'},'Location','best');
title('Verification Results per Configuration');
grid on;

% Figure 2: Solved percentage per config.
f2 = figure('Name','Solved Rate');
totalAll = numVerifAll + numFalsAll + numUnknownAll;
solvedPct = (numVerifAll + numFalsAll) ./ max(totalAll,1) * 100;
bar(solvedPct,'FaceColor',[0.2 0.5 0.8]);
set(gca,'XTickLabel',configLabels,'XTickLabelRotation',45,'FontSize',fontSize);
ylabel('Solved [%]');
ylim([0 100]);
title('Solved Rate per Configuration');
grid on;

% Figure 3: Total verification time per config.
f3 = figure('Name','Total Time');
bar(totalTimeAll,'FaceColor',[0.8 0.5 0.2]);
set(gca,'XTickLabel',configLabels,'XTickLabelRotation',45,'FontSize',fontSize);
ylabel('Total Time [s]');
title('Total Verification Time per Configuration');
grid on;

% Figure 4: Cumulative solved over time.
f4 = figure('Name','Cumulative Solved');
hold on;
colors = lines(numCombinations);
legendEntries = {};
for c = 1:numCombinations
    T = allResults{c};
    if isempty(T)
        continue;
    end
    isSolved = strcmp(T.results,'unsat') | strcmp(T.results,'sat');
    if ~any(isSolved)
        continue;
    end
    solvedTimes = sort(T.verifTimes(isSolved));
    cumSolved = (1:length(solvedTimes))';
    solvedTimes = [0; solvedTimes];
    cumSolved = [0; cumSolved];
    stairs(solvedTimes,cumSolved,'Color',colors(c,:),'LineWidth',1.5);
    legendEntries{end+1} = configLabels{c};
end
% Format axes and legend.
hold off;
set(gca,'FontSize',fontSize);
xlabel('Time [s]');
ylabel('Cumulative Solved Instances');
if ~isempty(legendEntries)
    legend(legendEntries,'Location','southeast');
end
title('Cumulative Solved Instances over Time');
grid on;

% Save .fig files and close figures to free memory.
saveas(f1,fullfile(plotspath,'verification_results.fig'));
saveas(f2,fullfile(plotspath,'solved_rate.fig'));
saveas(f3,fullfile(plotspath,'total_time.fig'));
saveas(f4,fullfile(plotspath,'cumulative_solved.fig'));
close([f1 f2 f3 f4]);
fprintf('Plots saved to: %s\n',plotspath);
catch plotErr
    fprintf('WARNING: Plot generation failed: %s\n',plotErr.message);
end

end

% ------------------------------ END OF CODE ------------------------------
