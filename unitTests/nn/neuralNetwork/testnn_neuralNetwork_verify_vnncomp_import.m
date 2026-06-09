function res = testnn_neuralNetwork_verify_vnncomp_import()
% testnn_neuralNetwork_verify_vnncomp_import - test that ONNX import and
%    VNNLib parsing succeeds for every supported benchmark configuration
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp_import()
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
% See also: neuralNetwork/verify, prepare_instance

% Authors:       Benedikt Kellner
% Written:       13-March-2026
% Last update:   26-March-2026 (dynamic year/benchmark discovery)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Supported benchNames (from prepare_instance.m).
supportedBenchNames = {'test','acasxu_2023','cifar100', ...
    'collins_rul_cnn_2023','collins_rul_cnn','cora','dist_shift_2023', ...
    'linearizenn','metaroom_2023','nn4sys_2023','safenlp','tinyimagenet', ...
    'tllverifybench_2023','mnist_fc','oval21','reach_prob_density', ...
    'rl_benchmarks','sri_resnet_a','sri_resnet_b'};

% Discover available VNN-COMP benchmark repositories.
dataDir = [CORAROOT '/examples/nn/vnncomp/data'];
yearDirs = dir(fullfile(dataDir, 'vnncomp*_benchmarks'));
yearDirs = yearDirs([yearDirs.isdir]);

if isempty(yearDirs)
    CORAwarning('CORA:nn','No VNN-COMP benchmark repositories found in %s. Skipping test.', dataDir);
    return;
end

% Save original directory and ensure we restore it.
origDir = pwd;
cleanupDir = onCleanup(@() cd(origDir));

% Change to the vnncomp directory (pipeline expects relative paths).
cd([CORAROOT '/examples/nn/vnncomp']);
addpath(pwd);

% Track whether we tested at least one benchmark.
numTested = 0;
unsupportedTested = false;

% =========================================================================
% Test supported benchmarks across all available years
% =========================================================================

for yi = 1:length(yearDirs)
    benchBase = fullfile(yearDirs(yi).folder, yearDirs(yi).name, 'benchmarks');
    if ~isfolder(benchBase)
        continue;
    end

    fprintf('\n=== %s ===\n', yearDirs(yi).name);

    % List all benchmark directories for this year.
    subdirs = dir(benchBase);
    subdirs = subdirs([subdirs.isdir]);
    subdirs = subdirs(~startsWith({subdirs.name}, '.'));

    for si = 1:length(subdirs)
        dirName = subdirs(si).name;
        benchName = aux_dirToBenchName(dirName, supportedBenchNames);

        % --- Test unsupported benchmark (first one found) ---
        if ~unsupportedTested && ~ismember(benchName, supportedBenchNames)
            unsupportedTested = aux_testUnsupported( ...
                benchBase, dirName, benchName);
            continue;
        end

        % Skip unsupported benchmarks.
        if ~ismember(benchName, supportedBenchNames)
            continue;
        end

        benchDir = fullfile(benchBase, dirName);
        instancesFile = fullfile(benchDir, 'instances.csv');
        if ~isfile(instancesFile)
            fprintf('Skipping %s: instances.csv not found.\n', dirName);
            continue;
        end

        % For nn4sys_2023, find first instance with a lindex model.
        if strcmp(benchName, 'nn4sys_2023')
            line = aux_findInstanceByPattern(instancesFile, 'lindex');
        else
            line = aux_readInstanceLine(instancesFile, 0);
        end
        if isempty(line)
            fprintf('Skipping %s: no suitable instance found.\n', dirName);
            continue;
        end

        % Parse CSV: onnxRelPath,vnnlibRelPath,timeout
        parts = strsplit(strtrim(line), ',');
        % Strip leading ./ from paths (some benchmarks use ./onnx/...).
        onnxRel = regexprep(strtrim(parts{1}), '^\.\/', '');
        vnnlibRel = regexprep(strtrim(parts{2}), '^\.\/', '');
        modelPath = fullfile(benchDir, onnxRel);
        vnnlibPath = fullfile(benchDir, vnnlibRel);

        fprintf('Testing import: %s (%s) ...\n', dirName, benchName);

        % 1. Run prepare_instance (suppress verbose output).
        [prepOutput, prepRes] = evalc('prepare_instance(benchName, modelPath, vnnlibPath)');
        assert(prepRes == 0, ...
            sprintf('prepare_instance failed for %s:\n%s', benchName, prepOutput));

        % 2. Check .mat file was created and is loadable.
        matFile = getInstanceFilename(benchName, modelPath, vnnlibPath);
        cleanupMat = onCleanup(@() aux_deleteIfExists(matFile));
        assert(isfile(matFile), ...
            sprintf('.mat file not created for %s.', benchName));

        data = load(matFile, 'nn', 'X0', 'specs', 'permuteDims');

        % 3. Evaluate center point.
        nn = data.nn;
        X0 = data.X0;
        xc = center(X0{1});
        if data.permuteDims
            inSize = nn.layers{1}.inputSize;
            % Check for collins_rul_cnn special case.
            [~,modelName,~] = getInstanceFilename(benchName,modelPath,vnnlibPath);
            if strcmp(benchName,'collins_rul_cnn_2023') ...
                    && ~strcmp(modelName{1},'NN_rul_full_window_40')
                permInSize = inSize;
            else
                permInSize = inSize([2 1 3]);
            end
            xc = reshape(permute(reshape(xc, permInSize), [2 1 3]), [], 1);
        end
        yc = nn.evaluate(xc);
        assert(~isempty(yc), ...
            sprintf('nn.evaluate returned empty for %s.', benchName));

        % 4. Check output dimension matches spec.
        specs = data.specs;
        spec1 = specs(1);
        if isa(spec1.set, 'halfspace')
            specDim = length(spec1.set.c);
        else
            specDim = size(spec1.set.A, 2);
        end
        assert(size(yc,1) == specDim, ...
            sprintf('Output dim %d does not match spec dim %d for %s.', ...
            size(yc,1), specDim, benchName));

        fprintf('  OK\n');
        numTested = numTested + 1;
    end
end

if numTested == 0
    CORAwarning('CORA:nn','No supported benchmarks found across any VNN-COMP year. Nothing tested.');
end

end


% Auxiliary functions -----------------------------------------------------

function benchName = aux_dirToBenchName(dirName, supportedNames)
    % Resolve directory name to a prepare_instance benchName by trying:
    %   1. dirName as-is
    %   2. dirName with year suffix stripped
    %   3. dirName with common year suffixes appended
    % Returns dirName unchanged if no supported match is found.

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

function tested = aux_testUnsupported(benchBase, dirName, benchName)
    % Test that prepare_instance returns 1 for an unsupported benchmark.
    tested = false;
    benchDir = fullfile(benchBase, dirName);
    instancesFile = fullfile(benchDir, 'instances.csv');
    if ~isfile(instancesFile)
        return;
    end
    line = aux_readInstanceLine(instancesFile, 0);
    if isempty(line)
        return;
    end
    parts = strsplit(strtrim(line), ',');
    onnxRel = regexprep(strtrim(parts{1}), '^\.\/', '');
    vnnlibRel = regexprep(strtrim(parts{2}), '^\.\/', '');
    modelPath = fullfile(benchDir, onnxRel);
    vnnlibPath = fullfile(benchDir, vnnlibRel);

    fprintf('Testing unsupported benchmark: %s (%s) ...\n', dirName, benchName);
    [~, prepRes] = evalc('prepare_instance(benchName, modelPath, vnnlibPath)');
    assert(prepRes == 1, ...
        'prepare_instance should return 1 for unsupported benchmark.');
    % Cleanup any .mat file that might have been created.
    matFile = getInstanceFilename(benchName, modelPath, vnnlibPath);
    aux_deleteIfExists(matFile);
    fprintf('  OK\n');
    tested = true;
end

function line = aux_readInstanceLine(instancesFile, lineNum)
    % Read a specific line from instances.csv.
    % lineNum = 0 means read the first line.
    fid = fopen(instancesFile, 'r');
    if fid == -1
        line = '';
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    if lineNum <= 0
        lineNum = 1;
    end
    line = '';
    for i = 1:lineNum
        line = fgetl(fid);
        if ~ischar(line)
            line = '';
            return;
        end
    end
end

function line = aux_findInstanceByPattern(instancesFile, pattern)
    % Find the first line in instances.csv containing the given pattern.
    fid = fopen(instancesFile, 'r');
    if fid == -1
        line = '';
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    line = '';
    while ~feof(fid)
        candidate = fgetl(fid);
        if ischar(candidate) && contains(candidate, pattern)
            line = candidate;
            return;
        end
    end
end

function aux_deleteIfExists(filepath)
    if isfile(filepath)
        delete(filepath);
    end
end

% ------------------------------ END OF CODE ------------------------------
