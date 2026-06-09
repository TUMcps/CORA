function res = testnn_neuralNetwork_verify_vnncomp()
% testnn_neuralNetwork_verify_vnncomp - end-to-end test of the vnncomp
%    pipeline (prepare_instance -> run_instance) with known ground truth
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp()
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
% See also: neuralNetwork/verify

% Authors:       Benedikt Kellner, Lukas Koller
% Written:       12-March-2026
% Last update:   26-March-2026 (dynamic year discovery)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Timeout for verification.
timeout = 30;
verbose = false;

% Discover available VNN-COMP benchmark repositories with a test benchmark.
dataDir = [CORAROOT '/examples/nn/vnncomp/data'];
yearDirs = dir(fullfile(dataDir, 'vnncomp*_benchmarks'));
yearDirs = yearDirs([yearDirs.isdir]);

testBenchDirs = {};
for i = 1:length(yearDirs)
    candidate = fullfile(yearDirs(i).folder, yearDirs(i).name, ...
        'benchmarks', 'test');
    if isfolder(candidate)
        testBenchDirs{end+1} = struct('path', candidate, ...
            'name', yearDirs(i).name);
    end
end

if isempty(testBenchDirs)
    CORAwarning('CORA:nn',['No VNN-COMP benchmark repositories with a test/ benchmark ' ...
        'found in %s. Skipping test.'], dataDir);
    return;
end

% Save original directory and ensure we restore it.
origDir = pwd;
cleanupDir = onCleanup(@() cd(origDir));

% Change to the vnncomp directory (pipeline expects relative paths).
cd([CORAROOT '/examples/nn/vnncomp']);
% Add vnncomp directory to path so prepare_instance/run_instance are found.
addpath(pwd);

benchName = 'test';

for ti = 1:length(testBenchDirs)
    testBenchDir = testBenchDirs{ti}.path;
    fprintf('\n=== %s ===\n', testBenchDirs{ti}.name);

    % Collect temp files for cleanup at end of this year.
    tempFiles = {};

    % =====================================================================
    % Test 1: test_sat (known SAT)
    % =====================================================================

    modelPath = [testBenchDir '/onnx/test_sat.onnx'];
    vnnlibPath = [testBenchDir '/vnnlib/test_prop.vnnlib'];

    if ~isfile(modelPath) || ~isfile(vnnlibPath)
        fprintf('Skipping test_sat: model or spec not found.\n');
    else
        resultsFile = [tempname '.txt'];
        tempFiles{end+1} = resultsFile;

        % Run prepare_instance.
        prepRes = prepare_instance(benchName,modelPath,vnnlibPath);
        assert(prepRes == 0, 'prepare_instance failed for test_sat.');

        % Check .mat file was created.
        matFile = getInstanceFilename(benchName,modelPath,vnnlibPath);
        assert(isfile(matFile), '.mat file not created by prepare_instance.');

        % Run run_instance (this deletes the .mat file internally).
        [resStr,~] = run_instance(benchName,modelPath,vnnlibPath, ...
            resultsFile,timeout,verbose);

        % Known SAT: must never claim unsat.
        assert(~strcmp(resStr,'unsat'), ...
            'test_sat: incorrectly returned unsat.');
        assert(strcmp(resStr,'sat') || strcmp(resStr,'unknown'), ...
            sprintf('test_sat: unexpected result ''%s''.',resStr));

        % Verify results file was created and content matches.
        assert(isfile(resultsFile), 'Results file not created.');
        fileContent = fileread(resultsFile);
        assert(contains(fileContent,resStr), ...
            'Results file content does not match return value.');

        % If sat, validate the counterexample from the results file.
        if strcmp(resStr,'sat')
            [x_parsed,y_parsed] = aux_parseCounterexample(resultsFile);
            % Reload network and spec.
            nn = neuralNetwork.readONNXNetwork(modelPath,false,'','', ...
                'dlnetwork',false);
            [X0, ~] = vnnlib2cora(vnnlibPath);
            % Check x_ within input bounds.
            tol = 1e-4; % %f format loses precision
            assert(all(x_parsed >= X0{1}.inf - tol) && ...
                all(x_parsed <= X0{1}.sup + tol), ...
                'Parsed counterexample x_ out of input bounds.');
            % Check output matches network evaluation.
            y_eval = nn.evaluate(x_parsed);
            assert(all(abs(y_eval - y_parsed) <= tol), ...
                'Parsed counterexample y_ does not match nn.evaluate(x_).');
        end
    end

    % =====================================================================
    % Test 2: test_unsat (known UNSAT)
    % =====================================================================

    modelPath2 = [testBenchDir '/onnx/test_unsat.onnx'];
    vnnlibPath2 = [testBenchDir '/vnnlib/test_prop.vnnlib'];

    if ~isfile(modelPath2) || ~isfile(vnnlibPath2)
        fprintf('Skipping test_unsat: model or spec not found.\n');
    else
        resultsFile2 = [tempname '.txt'];
        tempFiles{end+1} = resultsFile2;

        % Run prepare_instance.
        prepRes = prepare_instance(benchName,modelPath2,vnnlibPath2);
        assert(prepRes == 0, 'prepare_instance failed for test_unsat.');

        % Check .mat file was created.
        matFile2 = getInstanceFilename(benchName,modelPath2,vnnlibPath2);
        assert(isfile(matFile2), '.mat file not created for test_unsat.');

        % Run run_instance.
        [resStr2,~] = run_instance(benchName,modelPath2,vnnlibPath2, ...
            resultsFile2,timeout,verbose);

        % Known UNSAT: must never find a counterexample.
        assert(~strcmp(resStr2,'sat'), ...
            'test_unsat: incorrectly returned sat.');
        assert(strcmp(resStr2,'unsat') || strcmp(resStr2,'unknown'), ...
            sprintf('test_unsat: unexpected result ''%s''.',resStr2));

        % Verify results file.
        assert(isfile(resultsFile2), 'Results file not created for test_unsat.');
        fileContent2 = fileread(resultsFile2);
        assert(contains(fileContent2,resStr2), ...
            'Results file content does not match for test_unsat.');
    end

    % =====================================================================
    % Test 3: test_nano (simplest, smoke test)
    % =====================================================================

    modelPath3 = [testBenchDir '/onnx/test_nano.onnx'];
    vnnlibPath3 = [testBenchDir '/vnnlib/test_nano.vnnlib'];

    if ~isfile(modelPath3) || ~isfile(vnnlibPath3)
        fprintf('Skipping test_nano: model or spec not found.\n');
    else
        resultsFile3 = [tempname '.txt'];
        tempFiles{end+1} = resultsFile3;

        % Run prepare_instance.
        prepRes = prepare_instance(benchName,modelPath3,vnnlibPath3);
        assert(prepRes == 0, 'prepare_instance failed for test_nano.');

        matFile3 = getInstanceFilename(benchName,modelPath3,vnnlibPath3);
        assert(isfile(matFile3), '.mat file not created for test_nano.');

        % Run run_instance.
        [resStr3,~] = run_instance(benchName,modelPath3,vnnlibPath3, ...
            resultsFile3,timeout,verbose);

        % Smoke test: should return a valid result string.
        assert(ismember(resStr3,{'sat','unsat','unknown'}), ...
            sprintf('test_nano: invalid result ''%s''.',resStr3));

        % Verify results file exists.
        assert(isfile(resultsFile3), 'Results file not created for test_nano.');
    end

    % =====================================================================
    % Cleanup temp files
    % =====================================================================

    for i = 1:length(tempFiles)
        if isfile(tempFiles{i})
            delete(tempFiles{i});
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [x_,y_] = aux_parseCounterexample(resultsFile)
    % Parse counterexample values from a VNN-COMP results file.
    content = fileread(resultsFile);
    % Parse X values: (X_i value)
    xTokens = regexp(content,'\(X_\d+\s+([-\d.eE+]+)\)','tokens');
    x_ = cellfun(@(t) str2double(t{1}), xTokens)';
    % Parse Y values: (Y_i value)
    yTokens = regexp(content,'\(Y_\d+\s+([-\d.eE+]+)\)','tokens');
    y_ = cellfun(@(t) str2double(t{1}), yTokens)';
end

% ------------------------------ END OF CODE ------------------------------
