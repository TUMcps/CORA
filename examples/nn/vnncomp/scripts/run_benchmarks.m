function run_benchmarks(benchmarks,datapath,resultspath,options,instanceIds,timeoutMultiplier,preferredVersion)
% run_benchmarks - run all benchmarks.
%
% Syntax:
%    run_benchmarks(benchmarks,datapath,resultspath)
%    run_benchmarks(benchmarks,datapath,resultspath,options)
%    run_benchmarks(benchmarks,datapath,resultspath,options,instanceIds)
%    run_benchmarks(benchmarks,datapath,resultspath,options,instanceIds,timeoutMultiplier)
%    run_benchmarks(benchmarks,datapath,resultspath,options,instanceIds,timeoutMultiplier,preferredVersion)
%
% Inputs:
%    benchmarks - names of the benchmarks
%    datapath - path to the benchmark files
%    resultsPath - path to the results directory
%    options - (optional) fully resolved options struct. When omitted,
%              getDefaultVNNCOMPoptions is called per benchmark.
%    instanceIds - (optional) array of instance indices to run (default: all)
%    timeoutMultiplier - (optional) fraction of official timeout to use (default: 1.0)
%    preferredVersion - (optional) VNN-LIB version subdir to prefer when a
%              benchmark ships multiple (e.g. '1.0'/'2.0'); falls back to the
%              highest available version if the preferred one is absent
%              (default: '1.0')
%
% Outputs:
%    -
%
% References:
%    [1] VNN-COMP'24
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   06-June-2026 (BK, 2026 paths, versioned layout, preferredVersion flag)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Default: no explicit options (each benchmark loads its own defaults).
if nargin < 4
    options = [];
end
% Default: run all instances.
if nargin < 5
    instanceIds = [];
end
if nargin < 6 || isempty(timeoutMultiplier)
    timeoutMultiplier = 1.0;
end
% Default: prefer the 1.0 format (counterexamples are validated there).
if nargin < 7 || isempty(preferredVersion)
    preferredVersion = '1.0';
end

% Run all benchmarks in the current directory.

% Restrict number of CPU threads.
maxNumCompThreads(4);

% Get the basepath.
basepath = pwd;
% Make resultspath absolute (if not already).
if ~java.io.File(resultspath).isAbsolute()
    resultspath = sprintf('%s/%s',basepath,resultspath);
end
% Create a results directory.
mkdir(resultspath);

% Specify base directory.
benchpaths = {...
    sprintf('%s/vnncomp2026_benchmarks/benchmarks',datapath),...
    sprintf('%s/vnncomp2025_benchmarks/benchmarks',datapath),...
    sprintf('%s/vnncomp2024_benchmarks/benchmarks',datapath),...
    sprintf('%s/vnncomp2022_benchmarks/benchmarks',datapath),...
};

% List all benchmarks.
benchdirs = [];
for i=1:length(benchpaths)
    benchdirs = [benchdirs; dir(benchpaths{i})];
end
benchdirs = benchdirs( ...
    ~ismember({benchdirs.name}, {'.', '..'}) ...
    & ismember({benchdirs.name}, benchmarks));

% Warn about benchmarks that were not found.
if ~isempty(benchdirs)
    foundNames = {benchdirs.name};
else
    foundNames = {};
end
for i=1:length(benchmarks)
    if ~ismember(benchmarks{i}, foundNames)
        % Collect all available benchmark names across all data paths.
        allAvail = {};
        for j=1:length(benchpaths)
            d = dir(benchpaths{j});
            allAvail = [allAvail, setdiff({d.name},{'.','..','results'})];
        end
        fprintf('WARNING: Benchmark ''%s'' not found in data directories.\n', ...
            benchmarks{i});
        fprintf('  Searched: %s\n', strjoin(benchpaths, ', '));
        fprintf('  Available: %s\n', strjoin(unique(allAvail), ', '));
    end
end

for i=1:length(benchdirs)
    % Get the benchmark name.
    benchnamei = benchdirs(i).name;
    benchpathi = benchdirs(i).folder;
    fprintf('__________________________________________________________________\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('BENCHMARK %s (%d/%d)\n',benchnamei,i,length(benchdirs));
    fprintf('------------------------------------------------------------------\n');
    % Derive year from the data path (e.g. vnncomp2026_benchmarks -> '2026').
    m = regexp(benchpathi, 'vnncomp(\d{4})_benchmarks', 'tokens', 'once');
    if ~isempty(m)
        benchresultpath = sprintf('%s/%s_%s', resultspath, m{1}, benchnamei);
    else
        benchresultpath = sprintf('%s/%s', resultspath, benchnamei);
    end
    mkdir(benchresultpath);
    % Change directory to the current benchmark.
    cd(sprintf('%s/%s',benchpathi,benchnamei));
    % If the benchmark uses a versioned layout (e.g. 1.0/ and 2.0/ subdirs),
    % descend into one that contains instances.csv. Prefer preferredVersion
    % when present, otherwise fall back to the highest available version.
    % Use isfile() (filesystem-only) to avoid false positives from MATLAB path.
    if ~isfile('instances.csv')
        if isfile(fullfile(preferredVersion,'instances.csv'))
            cd(preferredVersion);
        else
            % collect versioned subdirs that contain instances.csv
            subdirs = dir('.');
            subdirs = subdirs([subdirs.isdir] & ~ismember({subdirs.name},{'.','..'}));
            subdirs = subdirs(arrayfun(@(s) ...
                isfile(fullfile(s.name,'instances.csv')), subdirs));
            if ~isempty(subdirs)
                % descend into the highest-versioned subdir (e.g. 2.0 over 1.0)
                [~,order] = sort({subdirs.name});
                cd(subdirs(order(end)).name);
            end
        end
    end
    % Resolve options for this benchmark: use provided options or load defaults.
    if isempty(options)
        benchOptions = getDefaultVNNCOMPoptions(benchnamei);
    else
        benchOptions = options;
    end
    % Run all instances of the benchmark.
    [numVerif,numFals,numUnknown] = ...
        run_instances(benchnamei,benchresultpath,benchOptions,instanceIds,timeoutMultiplier);
    % Compute total number of instances.
    totalNum = numVerif + numFals + numUnknown;
    % Print summary.
    table = CORAtableParameters(sprintf('BENCHMARK %s RESULTS',benchnamei));
    table.printHeader();
    table.printContentRow('#Verified',sprintf('%d/%d [%.1f%%]', ...
        numVerif,totalNum,numVerif/totalNum*100));
    table.printContentRow('#Falsified',sprintf('%d/%d [%.1f%%]', ...
        numFals,totalNum,numFals/totalNum*100));
    table.printContentRow('#Unknown',sprintf('%d/%d [%.1f%%]', ...
        numUnknown,totalNum,numUnknown/totalNum*100));
    table.printMidBoundaryRow();
    table.printContentRow('Solved',sprintf('%d/%d [%.1f%%]', ...
        numVerif+numFals,totalNum, ...
        (numVerif+numFals)/totalNum*100));
    table.printFooter();
    fprintf('__________________________________________________________________\n');
    % Go back to main directory.
    cd(basepath);
end

end

% ------------------------------ END OF CODE ------------------------------
