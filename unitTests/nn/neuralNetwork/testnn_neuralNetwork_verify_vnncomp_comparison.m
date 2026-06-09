function res = testnn_neuralNetwork_verify_vnncomp_comparison()
% testnn_neuralNetwork_verify_vnncomp_comparison - compare CORA verification
%    results against VNN-COMP competitor consensus ground truth
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp_comparison()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: prepare_instance.m, run_instance.m
% Subfunctions: aux_buildConsensus, aux_normalizeResult, aux_instanceKey
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Kellner
% Written:       12-March-2026
% Last update:   26-March-2026 (dynamic discovery, GitHub-fetched results)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Configuration -----------------------------------------------------------

numInstancesPerBench = 10;
minConsensus = 2; % >= 2 tools must agree for ground truth

% Supported benchNames for comparison (from prepare_instance.m, no 'test').
supportedBenchNames = {'acasxu_2023','cifar100','collins_rul_cnn_2023', ...
    'collins_rul_cnn','cora','dist_shift_2023','linearizenn', ...
    'metaroom_2023','nn4sys_2023','safenlp','tinyimagenet', ...
    'tllverifybench_2023','mnist_fc','oval21','reach_prob_density', ...
    'rl_benchmarks','sri_resnet_a','sri_resnet_b'};

% Competitor tools (exclude 'cora' - don't compare against ourselves).
competitors = {'alpha_beta_crown','neuralsat','pyrat','nnenum','nnv', ...
    'rover','marabou','sobolbox','never2','fastbatllnn'};

% Paths.
coraRoot = CORAROOT;
dataDir = fullfile(coraRoot,'examples','nn','vnncomp','data');
vnncompDir = fullfile(coraRoot,'examples','nn','vnncomp');

% Discover available VNN-COMP benchmark repositories ---------------------

yearDirs = dir(fullfile(dataDir, 'vnncomp*_benchmarks'));
yearDirs = yearDirs([yearDirs.isdir]);

if isempty(yearDirs)
    CORAwarning('CORA:nn',sprintf('No VNN-COMP benchmark repositories found in %s. Skipping test.', dataDir));
    res = true;
    return;
end

% Extract year numbers.
years = [];
for i = 1:length(yearDirs)
    tokens = regexp(yearDirs(i).name, 'vnncomp(\d{4})_benchmarks', 'tokens');
    if ~isempty(tokens)
        years(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
    end
end
years = sort(years);

% Build list of testable benchmarks: {year, dirName, benchName, benchPath}.
benchList = struct('year',{},'dirName',{},'benchName',{},'benchPath',{});
for yi = 1:length(years)
    yr = years(yi);
    benchBase = fullfile(dataDir, sprintf('vnncomp%d_benchmarks', yr), ...
        'benchmarks');
    if ~isfolder(benchBase)
        continue;
    end
    % Check subdirectories.
    subdirs = dir(benchBase);
    subdirs = subdirs([subdirs.isdir]);
    subdirs = subdirs(~startsWith({subdirs.name}, '.'));
    for si = 1:length(subdirs)
        dirName = subdirs(si).name;
        benchName = aux_dirToBenchName(dirName, supportedBenchNames);
        if ismember(benchName, supportedBenchNames) ...
                && isfile(fullfile(benchBase, dirName, 'instances.csv'))
            benchList(end+1) = struct('year', yr, 'dirName', dirName, ...
                'benchName', benchName, ...
                'benchPath', fullfile(benchBase, dirName)); %#ok<AGROW>
        end
    end
end

if isempty(benchList)
    CORAwarning('CORA:nn','No supported benchmarks found. Skipping test.');
    res = true;
    return;
end

% Store original directory and ensure we restore it.
origDir = pwd;
cleanupObj = onCleanup(@() cd(origDir));

% Add vnncomp directory to path so prepare_instance/run_instance are found.
addpath(vnncompDir);

% Step 1: Fetch consensus and run CORA per year ---------------------------

soundnessOk = true;
numCorrect = 0;
numMiss = 0;
numNoConsensus = 0;
numConflict = 0;
numSoundnessViolation = 0;

% Build dynamic benchmark stats using sanitized keys.
benchYears = arrayfun(@(b) b.year, benchList);
uniqueYears = unique(benchYears);
benchStatKeys = {};
benchStatLabels = {};
benchStats = struct();
for bi = 1:length(benchList)
    statKey = aux_statKey(benchList(bi).benchName, benchList(bi).year);
    if ~isfield(benchStats, statKey)
        benchStatKeys{end+1} = statKey; %#ok<AGROW>
        benchStatLabels{end+1} = sprintf('%s (%d)', ...
            benchList(bi).benchName, benchList(bi).year); %#ok<AGROW>
        benchStats.(statKey) = struct( ...
            'correct',0,'miss',0,'noConsensus',0,'conflict',0, ...
            'violation',0,'total',0);
    end
end

allResults = struct('statKey',{},'onnxRel',{},'vnnlibRel',{}, ...
    'coraResult',{});

for yi = 1:length(uniqueYears)
    yr = uniqueYears(yi);

    % --- Fetch competitor consensus from GitHub ---
    fprintf('\n=== Year %d: Fetching competitor consensus from GitHub ===\n', yr);
    consensus = aux_buildConsensus(yr, competitors, minConsensus);

    if consensus.Count == 0
        fprintf('No consensus data for year %d. Skipping.\n', yr);
        continue;
    end
    fprintf('Consensus computed for %d instances.\n', consensus.Count);

    % --- Run CORA on benchmarks for this year ---
    yearBenches = benchList(benchYears == yr);

    for bi = 1:length(yearBenches)
        bench = yearBenches(bi);
        statKey = aux_statKey(bench.benchName, bench.year);

        fprintf('\n--- Benchmark: %s (%d) ---\n', bench.benchName, bench.year);

        % Read instances.csv.
        instancesFile = fullfile(bench.benchPath, 'instances.csv');
        fid = fopen(instancesFile, 'r');
        if fid == -1
            CORAwarning('CORA:nn',sprintf('Cannot open instances.csv for "%s". Skipping.', ...
                bench.benchName));
            continue;
        end
        rawLines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        rawLines = rawLines{1};

        numInstances = min(numInstancesPerBench, length(rawLines));

        for i = 1:numInstances
            parts = strsplit(rawLines{i}, ',');
            modelRel = strtrim(parts{1});
            vnnlibRel = strtrim(parts{2});
            instanceTimeout = str2double(strtrim(parts{3}));

            fprintf('\n  Instance %d/%d: %s | %s (timeout=%ds)\n', ...
                i, numInstances, modelRel, vnnlibRel, instanceTimeout);

            % Change to benchmark directory (prepare_instance uses
            % relative paths).
            cd(bench.benchPath);

            % Results file for this instance.
            [~,modelBase] = fileparts(modelRel);
            [~,vnnlibBase] = fileparts(vnnlibRel);
            resultsFile = fullfile(bench.benchPath, ...
                sprintf('result_%s_%s.txt', modelBase, vnnlibBase));

            coraResult = 'unknown';
            try
                prepRes = prepare_instance(bench.benchName, ...
                    modelRel, vnnlibRel);
                if prepRes ~= 0
                    fprintf('  prepare_instance failed (code=%d). Treating as unknown.\n', prepRes);
                else
                    [resStr,~] = run_instance(bench.benchName, ...
                        modelRel, vnnlibRel, resultsFile, ...
                        instanceTimeout, false);
                    coraResult = resStr;
                end
            catch e
                fprintf('  Error: %s\n', e.message);
                coraResult = 'unknown';
            end

            % Clean up temp files.
            aux_cleanupTempFiles(bench.benchName, modelRel, ...
                vnnlibRel, resultsFile);

            fprintf('  CORA result: %s\n', coraResult);

            entry.statKey = statKey;
            entry.onnxRel = modelRel;
            entry.vnnlibRel = vnnlibRel;
            entry.coraResult = coraResult;
            allResults(end+1) = entry; %#ok<AGROW>
        end
    end

    % --- Compare results for this year ---
    fprintf('\n=== Year %d: Comparing results ===\n', yr);

    for k = 1:length(allResults)
        r = allResults(k);
        % Only process results from this year.
        if ~startsWith(r.statKey, sprintf('y%d_', yr))
            continue;
        end

        key = aux_instanceKey(aux_extractRelPath(r.onnxRel), ...
            aux_extractRelPath(r.vnnlibRel));
        sk = r.statKey;
        benchStats.(sk).total = benchStats.(sk).total + 1;

        if ~consensus.isKey(key)
            numNoConsensus = numNoConsensus + 1;
            benchStats.(sk).noConsensus = benchStats.(sk).noConsensus + 1;
            fprintf('  [no consensus] %s: CORA=%s\n', key, r.coraResult);
            continue;
        end

        cons = consensus(key);

        if strcmp(cons.status, 'conflict')
            numConflict = numConflict + 1;
            benchStats.(sk).conflict = benchStats.(sk).conflict + 1;
            fprintf('  [conflict] %s: CORA=%s (sat=%d, unsat=%d)\n', ...
                key, r.coraResult, cons.numSat, cons.numUnsat);
            continue;
        end

        gt = cons.result;

        if strcmp(r.coraResult, 'unknown')
            numMiss = numMiss + 1;
            benchStats.(sk).miss = benchStats.(sk).miss + 1;
            fprintf('  [miss] %s: CORA=unknown, consensus=%s\n', key, gt);
        elseif strcmp(r.coraResult, gt)
            numCorrect = numCorrect + 1;
            benchStats.(sk).correct = benchStats.(sk).correct + 1;
            fprintf('  [correct] %s: CORA=%s, consensus=%s\n', ...
                key, r.coraResult, gt);
        else
            numSoundnessViolation = numSoundnessViolation + 1;
            benchStats.(sk).violation = benchStats.(sk).violation + 1;
            soundnessOk = false;
            fprintf('  [VIOLATION] %s: CORA=%s, consensus=%s\n', ...
                key, r.coraResult, gt);
        end
    end
end

% Restore directory.
cd(origDir);

% Step 2: Report ----------------------------------------------------------

fprintf('\n=== Summary ===\n\n');

fprintf('%-35s %8s %8s %8s %8s %8s %8s\n', ...
    'Benchmark','Total','Correct','Miss','NoCons','Conflict','Violat.');
fprintf('%s\n', repmat('-',1,93));
for bi = 1:length(benchStatKeys)
    sk = benchStatKeys{bi};
    s = benchStats.(sk);
    fprintf('%-35s %8d %8d %8d %8d %8d %8d\n', ...
        benchStatLabels{bi}, s.total, s.correct, s.miss, ...
        s.noConsensus, s.conflict, s.violation);
end
fprintf('%s\n', repmat('-',1,93));
fprintf('%-35s %8d %8d %8d %8d %8d %8d\n', ...
    'TOTAL', length(allResults), numCorrect, numMiss, numNoConsensus, ...
    numConflict, numSoundnessViolation);

fprintf('\nLegend:\n');
fprintf('  Correct  - CORA result matches competitor consensus\n');
fprintf('  Miss     - CORA returned unknown, but consensus exists (timeout/error)\n');
fprintf('  NoCons   - Not enough competitor votes to establish ground truth\n');
fprintf('  Conflict - Competitors disagree (some say sat, others unsat)\n');
fprintf('  Violat.  - CORA contradicts consensus (soundness violation!)\n');
fprintf('\nSoundness OK: %d\n', soundnessOk);

% Assert no soundness violations.
assert(soundnessOk, ...
    'Soundness violation detected! CORA contradicts competitor consensus.');

res = true;

end


% Auxiliary functions -----------------------------------------------------

function benchName = aux_dirToBenchName(dirName, supportedNames)
% aux_dirToBenchName - resolve directory name to prepare_instance benchName
%    by trying: exact match, year suffix stripped, year suffixes appended.
    % 1. Exact match.
    if ismember(dirName, supportedNames)
        benchName = dirName;
        return;
    end
    % 2. Strip trailing _YYYY suffix.
    base = regexprep(dirName, '_\d{4}$', '');
    if ~strcmp(base, dirName) && ismember(base, supportedNames)
        benchName = base;
        return;
    end
    % 3. Try appending year suffixes to the (possibly stripped) base.
    yearSuffixes = {'_2023','_2024','_2025','_2022'};
    for i = 1:length(yearSuffixes)
        candidate = [base yearSuffixes{i}];
        if ismember(candidate, supportedNames)
            benchName = candidate;
            return;
        end
    end
    % 4. Fallback: return dirName unchanged.
    benchName = dirName;
end


function sk = aux_statKey(benchName, year)
% aux_statKey - create a valid MATLAB struct fieldname from benchmark+year.
    sk = sprintf('y%d_%s', year, regexprep(benchName, '\W', '_'));
end


function consensus = aux_buildConsensus(year, tools, minConsensus)
% aux_buildConsensus - fetch competitor CSVs from GitHub and compute
%    per-instance consensus ground truth.
%
%    Returns a containers.Map from instance key to a struct with fields:
%      .result   - 'sat' or 'unsat' (only if consensus reached)
%      .status   - 'consensus' or 'conflict'
%      .numSat   - number of tools that said sat
%      .numUnsat - number of tools that said unsat

    votes = containers.Map('KeyType','char','ValueType','any');
    % Specify VNNCOMP results url.
    urlBase = sprintf( ...
        'https://raw.githubusercontent.com/VNN-COMP/vnncomp%d_results/main', ...
        year);
    opts = weboptions('ContentType','text','Timeout',30);

    % Check the results of all tools.
    for i = 1:length(tools)
        toolName = tools{i};
        url = sprintf('%s/%s/results.csv', urlBase, toolName);

        % Read results file.
        try
            csvText = webread(url, opts);
        catch
            fprintf('  [fetch failed] %s - skipping\n', toolName);
            continue;
        end

        fprintf('  Fetched %s ...', toolName);
        rawLines = strsplit(csvText, newline);

        % Check all lines.
        numEntries = 0;
        for j = 1:length(rawLines)
            line = rawLines{j};
            if isempty(strtrim(line))
                continue;
            end
            parts = strsplit(line, ',');
            if length(parts) < 5
                continue;
            end

            % Extract the model, vnnlib, and result.
            onnxPath = strtrim(parts{2});
            vnnlibPath = strtrim(parts{3});
            resultStr = strtrim(parts{5});

            normResult = aux_normalizeResult(resultStr);
            if strcmp(normResult, 'inconclusive')
                continue;
            end

            onnxRel = aux_extractRelPath(onnxPath);
            vnnlibRel = aux_extractRelPath(vnnlibPath);
            key = aux_instanceKey(onnxRel, vnnlibRel);

            if votes.isKey(key)
                v = votes(key);
            else
                v = struct('numSat',0,'numUnsat',0);
            end

            % Aggregate result.
            if strcmp(normResult,'sat')
                v.numSat = v.numSat + 1;
            elseif strcmp(normResult,'unsat')
                v.numUnsat = v.numUnsat + 1;
            end

            votes(key) = v;
            numEntries = numEntries + 1;
        end
        fprintf(' %d entries\n', numEntries);
    end

    % Build consensus from votes.
    consensus = containers.Map('KeyType','char','ValueType','any');
    keys = votes.keys();
    for i = 1:length(keys)
        key = keys{i};
        v = votes(key);

        % Construct entry.
        entry = struct('result','','status','', ...
            'numSat',v.numSat,'numUnsat',v.numUnsat);

        % Set entry data.
        if v.numSat > 0 && v.numUnsat > 0
            entry.status = 'conflict';
            entry.result = '';
        elseif v.numSat >= minConsensus && v.numUnsat == 0
            entry.status = 'consensus';
            entry.result = 'sat';
        elseif v.numUnsat >= minConsensus && v.numSat == 0
            entry.status = 'consensus';
            entry.result = 'unsat';
        else
            continue;
        end

        consensus(key) = entry;
    end
end


function normResult = aux_normalizeResult(result)
% aux_normalizeResult - map VNN-COMP result strings to sat/unsat/inconclusive.
    result = lower(strtrim(result));
    if strcmp(result,'sat') || strcmp(result,'violated')
        normResult = 'sat';
    elseif strcmp(result,'unsat') || strcmp(result,'holds') || strcmp(result,'verified')
        normResult = 'unsat';
    else
        normResult = 'inconclusive';
    end
end


function key = aux_instanceKey(onnxRel, vnnlibRel)
% aux_instanceKey - create a unique key from onnx and vnnlib relative paths.
    key = [onnxRel '|' vnnlibRel];
end


function relPath = aux_extractRelPath(fullPath)
% aux_extractRelPath - extract the relative path starting from onnx/ or
%    vnnlib/ directory component.
    fullPath = strrep(fullPath, '\', '/');
    if startsWith(fullPath, './')
        fullPath = fullPath(3:end);
    end
    idx = regexp(fullPath, '(onnx/|vnnlib/)', 'once');
    if ~isempty(idx)
        relPath = fullPath(idx:end);
    else
        relPath = fullPath;
    end
end


function aux_cleanupTempFiles(benchName, modelRel, vnnlibRel, resultsFile)
% aux_cleanupTempFiles - remove temporary .mat and results files.
    try
        matFile = getInstanceFilename(benchName, modelRel, vnnlibRel);
        if isfile(matFile)
            delete(matFile);
        end
    catch
    end
    try
        if isfile(resultsFile)
            delete(resultsFile);
        end
    catch
    end
end

% ------------------------------ END OF CODE ------------------------------
