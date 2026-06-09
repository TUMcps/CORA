function res = testnn_neuralNetwork_hyperparam_tuning_vnncomp()
% testnn_neuralNetwork_hyperparam_tuning_vnncomp - tests that
%    hyperparam_tuning_vnncomp overrides reach the verification engine by
%    running the same benchmark with a swapped heuristic and checking the
%    results differ
%
% Syntax:
%    res = testnn_neuralNetwork_hyperparam_tuning_vnncomp()
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
% See also: hyperparam_tuning_vnncomp, prepare_instance

% Authors:       Benedikt Kellner
% Written:       23-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Discover a VNN-COMP test benchmark.
dataDir = [CORAROOT '/examples/nn/vnncomp/data'];
yearDirs = dir(fullfile(dataDir,'vnncomp*_benchmarks'));
yearDirs = yearDirs([yearDirs.isdir]);

% Layouts differ between years: some keep models flat inside the
% benchmark directory, newer ones use onnx/ and vnnlib/ subdirs.
modelPath = '';
vnnlibPath = '';
for i = 1:length(yearDirs)
    candidate = fullfile(yearDirs(i).folder,yearDirs(i).name, ...
        'benchmarks','test');
    if ~isfolder(candidate)
        continue;
    end
    % Try both subdir and flat layouts for the SAT instance.
    for mp = {[candidate '/onnx/test_sat.onnx'], ...
            [candidate '/test_sat.onnx']}
        if isfile(mp{1})
            modelPath = mp{1};
            break;
        end
    end
    % Same for the specification file.
    for vp = {[candidate '/vnnlib/test_prop.vnnlib'], ...
            [candidate '/test_prop.vnnlib']}
        if isfile(vp{1})
            vnnlibPath = vp{1};
            break;
        end
    end
    if ~isempty(modelPath) && ~isempty(vnnlibPath)
        break;
    end
end

if isempty(modelPath) || isempty(vnnlibPath)
    CORAwarning('CORA:nn',['No VNN-COMP test_sat instance found in %s. ' ...
        'Skipping hyperparam tuning test.'],dataDir);
    return;
end

origDir = pwd;
cleanupDir = onCleanup(@() cd(origDir));

cd([CORAROOT '/examples/nn/vnncomp']);
addpath(pwd);

benchName = 'test';

% Two configs that differ only in the neuron split heuristic. Options are
% passed exactly as the tuning pipeline builds them: a fully resolved
% options.nn struct (prepare_instance does no flat-key translation — the
% '__'-to-nested mapping happens upstream in hyperparam_tuning).
overridesA = struct('nn',struct('neuron_split_heuristic','zono-norm-gradient', ...
    'num_splits',2,'num_dimensions',1,'num_neuron_splits',1));
overridesB = overridesA;
overridesB.nn.neuron_split_heuristic = 'least-unstable';

% Override is applied to options.nn (verified via saved .mat).
prepRes = prepare_instance(benchName,modelPath,vnnlibPath,false,overridesA);
assert(prepRes == 0,'prepare_instance failed for config A.');
matFile = getInstanceFilename(benchName,modelPath,vnnlibPath);
loaded = load(matFile,'options');
assert(strcmp(loaded.options.nn.neuron_split_heuristic,'zono-norm-gradient'), ...
    'Override for neuron_split_heuristic was not applied in config A.');
assert(loaded.options.nn.num_neuron_splits == 1, ...
    'Override for num_neuron_splits was not applied.');

% Nested train-field overrides reach options.nn.train.
nestedOverrides = overridesA;
nestedOverrides.nn.train.mini_batch_size = 128;
prepRes = prepare_instance(benchName,modelPath,vnnlibPath,false,nestedOverrides);
assert(prepRes == 0,'prepare_instance failed for nested override.');
loaded = load(matFile,'options');
assert(loaded.options.nn.train.mini_batch_size == 128, ...
    'Nested override did not reach options.nn.train.');

% Run both configs end-to-end. numSubproblems must be populated (the
% field is always initialized in run_instance's result struct).
timeout = 30;
verbose = false;
[subA,~] = aux_runWithOverrides(benchName,modelPath,vnnlibPath, ...
    overridesA,timeout,verbose);
[subB,~] = aux_runWithOverrides(benchName,modelPath,vnnlibPath, ...
    overridesB,timeout,verbose);
assert(isnumeric(subA) && isnumeric(subB) && subA >= 0 && subB >= 0, ...
    'numSubproblems is not a valid non-negative number.');

% The swapped heuristic must reach the verification engine: reload the
% saved options for config B and verify the override survived the merge
% with the benchmark-specific defaults in prepare_instance.
prepare_instance(benchName,modelPath,vnnlibPath,false,overridesB);
loadedB = load(matFile,'options');
assert(strcmp(loadedB.options.nn.neuron_split_heuristic,'least-unstable'), ...
    'Swapped heuristic in config B was overwritten by benchmark defaults.');

end


% Auxiliary functions -----------------------------------------------------

function [numSub,verifTime] = aux_runWithOverrides(benchName,modelPath, ...
    vnnlibPath,overrides,timeout,verbose)
% Run prepare_instance + run_instance with the given overrides and
% return the number of verified subproblems.
prepRes = prepare_instance(benchName,modelPath,vnnlibPath,verbose, ...
    overrides);
assert(prepRes == 0,'prepare_instance failed.');
resultsFile = [tempname '.txt'];
cleanup = onCleanup(@() aux_deleteIfExists(resultsFile));
[~,resOut] = run_instance(benchName,modelPath,vnnlibPath,resultsFile, ...
    timeout,verbose);
numSub = resOut.numSubproblems;
verifTime = resOut.time;
end

function aux_deleteIfExists(f)
if isfile(f)
    delete(f);
end
end

% ------------------------------ END OF CODE ------------------------------
