function res = testnn_neuralNetwork_verify_vnncomp_multinet()
% testnn_neuralNetwork_verify_vnncomp_multinet - end-to-end test of the
%    multi-network pipeline (prepare_instance -> run_instance) on a
%    monotonic_acasxu_2026 instance: builds the joint network, falsifies the
%    (false) monotonicity property ('sat'), and checks the counterexample
%    re-evaluates consistently. Skips if the benchmark data is absent.
%
% Syntax:
%    res = testnn_neuralNetwork_verify_vnncomp_multinet()
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
% See also: neuralNetwork/verify, prepare_instance, run_instance,
%           nnHelper.buildJointNetwork

% Authors:       Benedikt Kellner
% Written:       03-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

verbose = false;
timeout = 60;
benchName = 'monotonic_acasxu_2026';

% Locate the benchmark directory (the one containing instances.csv).
benchRoot = [CORAROOT '/examples/nn/vnncomp/data/vnncomp2026_benchmarks' ...
    '/benchmarks/' benchName];
benchDir = aux_findInstancesDir(benchRoot);
if isempty(benchDir)
    CORAwarning('CORA:nn', ['monotonic_acasxu_2026 benchmark not found ' ...
        'under %s. Skipping multi-net end-to-end test.'], benchRoot);
    return;
end

% Save and restore working directory + path.
origDir = pwd;
cleanupDir = onCleanup(@() cd(origDir));
addpath([CORAROOT '/examples/nn/vnncomp']);
% Pipeline expects relative onnx/vnnlib paths -> run from the benchmark dir.
cd(benchDir);

% Read the first instance from instances.csv.
[modelPath, vnnlibPath] = aux_firstInstance(fullfile(benchDir,'instances.csv'));
if isempty(modelPath)
    CORAwarning('CORA:nn','instances.csv empty/unreadable. Skipping.');
    return;
end

% =====================================================================
% 1. prepare_instance builds the joint network + reduced input polytope
% =====================================================================
prepRes = prepare_instance(benchName, modelPath, vnnlibPath, verbose);
assert(prepRes == 0, 'prepare_instance failed for monotonic instance.');

matFile = getInstanceFilename(benchName, modelPath, vnnlibPath);
assert(isfile(matFile), '.mat file not created by prepare_instance.');

S = load(matFile);
% Joint network: a single concat composite layer.
assert(isa(S.nn, 'neuralNetwork'), 'stored nn is not a neuralNetwork.');
assert(isa(S.nn.layers{1}, 'nnCompositeLayer'), ...
    'joint network first layer is not an nnCompositeLayer.');
% Multi-network metadata.
assert(~isempty(S.multiNetInfo), 'multiNetInfo not stored.');
mni = S.multiNetInfo;
assert(isequal(mni.netNames, {'f','g'}), 'unexpected network names.');
assert(mni.nfOut == 5 && mni.ngOut == 5, 'ACAS Xu outputs should be 5+5.');
% Reduced input set: 6 free variables (5 X_f + 1 uncoupled X_g[0]).
assert(isscalar(S.X0), 'X0 should be a single set.');
assert(dim(S.X0{1}) == mni.nFree, 'X0 dim != nFree.');
assert(mni.nFree == 6, 'monotonic free-variable count should be 6.');
assert(isa(S.X0{1}, 'polytope'), ...
    'monotonic X0 should be a polytope (residual X_g0 <= X_f0 coupling).');
assert(isequal(size(mni.Sf), [5 6]) && isequal(size(mni.Sg), [5 6]), ...
    'selector matrices have wrong shape.');

% =====================================================================
% 2. run_instance falsifies the (false) monotonicity property
% =====================================================================
resultsFile = [tempname '.txt'];
cleanupRes = onCleanup(@() aux_safeDelete(resultsFile));

[resStr,~] = run_instance(benchName, modelPath, vnnlibPath, ...
    resultsFile, timeout, verbose);

% The property is genuinely false for ACAS Xu -> must never claim 'unsat'.
assert(~strcmp(resStr,'unsat'), ...
    'monotonic instance incorrectly verified (unsat).');
assert(strcmp(resStr,'sat'), ...
    sprintf('monotonic instance expected ''sat'', got ''%s''.', resStr));
assert(isfile(resultsFile), 'results file not created.');

% =====================================================================
% 3. Multi-network counterexample is consistent with the joint network
% =====================================================================
[xf, xg, yf, yg] = aux_parseMultiNetCE(resultsFile);
assert(numel(xf) == 5 && numel(xg) == 5, 'CE input split has wrong size.');
assert(numel(yf) == 5 && numel(yg) == 5, 'CE output split has wrong size.');

% Reconstruct the free-variable vector v:  v(1:nf)=xf,  v(nf+k)=xg(uncoupledG(k)).
nf = size(mni.Sf, 1);
v = zeros(mni.nFree, 1);
v(1:nf) = xf;
for k = 1:numel(mni.uncoupledG)
    v(nf + k) = xg(mni.uncoupledG(k));
end
% Selector matrices must reproduce the per-network inputs from v.
assert(all(withinTol(mni.Sf * v, xf(:), 1e-4)), 'Sf*v != xf.');
assert(all(withinTol(mni.Sg * v, xg(:), 1e-4)), 'Sg*v != xg.');
% Joint network evaluation must reproduce the reported [yf; yg].
y = S.nn.evaluate(v);
assert(all(withinTol(y, [yf(:); yg(:)], 1e-3)), ...
    'joint network output does not match the counterexample outputs.');

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function d = aux_findInstancesDir(benchRoot)
% Return the directory under benchRoot that contains instances.csv, or '' .
d = '';
if ~isfolder(benchRoot), return; end
if isfile(fullfile(benchRoot,'instances.csv'))
    d = benchRoot; return;
end
sub = dir(benchRoot);
for k = 1:numel(sub)
    if sub(k).isdir && ~startsWith(sub(k).name,'.') ...
            && isfile(fullfile(benchRoot, sub(k).name, 'instances.csv'))
        d = fullfile(benchRoot, sub(k).name); return;
    end
end
end

function [modelPath, vnnlibPath] = aux_firstInstance(csvFile)
% Parse the first line of an instances.csv:
%   "[('f','...'),('g','...')]",./vnnlib/instance_0.vnnlib,100
modelPath = ''; vnnlibPath = '';
if ~isfile(csvFile), return; end
fid = fopen(csvFile,'r');
line = fgetl(fid);
fclose(fid);
if ~ischar(line), return; end
tok = regexp(line, '^"(.*)",([^,]+),', 'tokens', 'once');
if isempty(tok), return; end
modelPath  = tok{1};
vnnlibPath = tok{2};
end

function [xf, xg, yf, yg] = aux_parseMultiNetCE(resultsFile)
% Parse a multi-network VNN-LIB 2.0 counterexample written by run_instance:
%   sat
%   <name> <dtype> [<shape>]
%   <value>   (one per line)
%   ...
content = fileread(resultsFile);
xf = aux_grab(content, 'X_f');
xg = aux_grab(content, 'X_g');
yf = aux_grab(content, 'Y_f');
yg = aux_grab(content, 'Y_g');
end

function vals = aux_grab(content, name)
% Collect the values in the block headed by "<name> <dtype> [<shape>]", up to
% the next header line or end of file.
lines = strsplit(strtrim(content), newline);
vals = [];
collecting = false;
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line), continue; end
    hdr = regexp(line, '^(\S+)\s+\w+\s+\[[\d,]*\]$', 'tokens', 'once');
    if ~isempty(hdr)
        collecting = strcmp(hdr{1}, name);
        continue;
    end
    if collecting
        v = str2double(line);
        if ~isnan(v), vals(end+1) = v; end %#ok<AGROW>
    end
end
vals = vals(:);
end

function aux_safeDelete(f)
if isfile(f), delete(f); end
end

% ------------------------------ END OF CODE ------------------------------
