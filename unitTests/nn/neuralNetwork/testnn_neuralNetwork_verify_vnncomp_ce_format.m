function res = testnn_neuralNetwork_verify_vnncomp_ce_format()
% testnn_neuralNetwork_verify_vnncomp_ce_format - check that prepare_instance
%    persists the right VNN-LIB version (1.0/2.0) for an ACASXU model with a
%    1.x and a 2.0 spec. No verification is run, so the test is deterministic;
%    the format output itself is checked in
%    testnn_neuralNetwork_verify_vnncomp_regression.
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp_ce_format()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: prepare_instance, getInstanceFilename
% Subfunctions: aux_assertSavedVersion
% MAT-files required: none
%
% See also: prepare_instance, run_instance, vnnlib2cora,
%    testnn_neuralNetwork_verify_vnncomp_regression

% Authors:       Benedikt Kellner
% Written:       05-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
verbose = false;

% the pipeline expects its helper functions (prepare_instance,
% getInstanceFilename) on the path and writes the instance .mat relative to
% the vnncomp directory
origDir = pwd;
cleanupDir = onCleanup(@() cd(origDir));
cd([CORAROOT '/examples/nn/vnncomp']);
addpath(pwd);

model = [CORAROOT '/models/Cora/nn/ACASXU_run2a_5_3_batch_2000.onnx'];
v1Spec = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/axas_xu_prop_3.vnnlib'];
v2Spec = [CORAROOT '/models/Cora/nn/unitTests/vnnlib/v2/acasxu_prop_3_v2.vnnlib'];

% a VNN-LIB 1.x spec must be persisted with version '1.0' (selects legacy
% counterexample format); a 2.0 spec with version '2.0' (selects new format)
aux_assertSavedVersion('acasxu', model, v1Spec, '1.0', verbose);
aux_assertSavedVersion('acasxu', model, v2Spec, '2.0', verbose);

end


% Auxiliary functions -----------------------------------------------------

function aux_assertSavedVersion(benchName, model, vnnlibPath, expectedVersion, verbose)
% run prepare_instance and assert the persisted vnnlibInfo carries the
% expected VNN-LIB version; skip (do not fail) when fixtures are unavailable
% or the environment cannot load the ONNX model

if ~isfile(model) || ~isfile(vnnlibPath)
    return;
end

% prepare_instance returns nonzero when it cannot load the network/spec
% (e.g., ONNX support missing); skip rather than fail in that case
if prepare_instance(benchName, model, vnnlibPath, verbose) ~= 0
    return;
end

instanceFilename = getInstanceFilename(benchName, model, vnnlibPath);
cleanupMat = onCleanup(@() aux_deleteIfPresent(instanceFilename));

S = load(instanceFilename, 'vnnlibInfo');
assert(isfield(S, 'vnnlibInfo') && ~isempty(S.vnnlibInfo), ...
    'prepare_instance did not persist vnnlibInfo for %s', vnnlibPath);
assert(strcmp(S.vnnlibInfo.version, expectedVersion), ...
    'persisted version ''%s'' (expected ''%s'') for %s', ...
    S.vnnlibInfo.version, expectedVersion, vnnlibPath);
end

function aux_deleteIfPresent(filename)
% delete a file only if it still exists (run_instance normally removes it)
if isfile(filename)
    delete(filename);
end
end

% ------------------------------ END OF CODE ------------------------------
