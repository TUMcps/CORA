function res = prepare_instance(benchName,modelPath,vnnlibPath,varargin)
% prepare_instance - called before the verification timer is started. Load
% the neural network and the vnnlib specification; they are stored as a
% .mat file, which is opened by run_instance.m
%
% Syntax:
%    res = prepare_instance(benchName,modelPath,vnnlibPath)
%    res = prepare_instance(benchName,modelPath,vnnlibPath,verbose)
%    res = prepare_instance(benchName,modelPath,vnnlibPath,verbose,options)
%
% Inputs:
%    benchName - name of the benchmark
%    modelPath - path to the .onnx-file
%    vnnlibPath - path to the .vnnlib-file
%    verbose - (optional) print progress, default true
%    options - (optional) fully resolved options struct. When omitted,
%              getDefaultVNNCOMPoptions(benchName) is used (competition / shell
%              script path). For tuning, always pass explicit options.
%
% Outputs:
%    res - result code
%
% References:
%    [1] VNN-COMP'24
%    [2] VNN-COMP'25
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getDefaultVNNCOMPoptions, run_instances

% Authors:       Lukas Koller, Benedikt Kellner
% Written:       11-August-2025
% Last update:   04-May-2026
%                06-June-2026 (BK, multi-network support and v1/v2 counterexample format)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[verbose,options] = setDefaultValues({true,[]},varargin);
if isempty(options)
    options = getDefaultVNNCOMPoptions(benchName);
end

% Fill in any missing required fields without overwriting paramGrid values.
options = nnHelper.validateNNoptions(options,true);
% Always override GPU detection — not a hyperparameter.
options.nn.train.use_gpu = aux_isGPUavailable();

% Parse the benchmark name.
[benchName_,~] = aux_parseBenchmarkName(benchName);

if verbose
    fprintf('prepare_instance(%s,%s,%s)...\n',benchName,modelPath, ...
        vnnlibPath);
end

try
    % The following steps are independent of the network type.
    % Resolve .gz compressed vnnlib if the plain path does not exist.
    vnnlibPath = aux_resolveVnncompPath(vnnlibPath);

    if verbose
        fprintf('--- GPU available: %d\n',options.nn.train.use_gpu);
        fprintf('--- Loading specification...');
    end
    % vnnlibInfo (version + per-network metadata) drives counterexample
    % formatting and, for multi-network specs, the joint-network build.
    [X0,specs,vnnlibInfo] = vnnlib2cora(vnnlibPath);
    if verbose
        fprintf(' done\n');
    end

    if ismember(benchName_,{'monotonic_acasxu','isomorphic_acasxu'})
        % ---- Multi-network path ------------------------------------------
        % Build the composite network from the parsed joint spec: the joint
        % polytope X0 encodes which input dims are coupled across networks.
        if verbose
            fprintf('--- Building joint network...');
        end
        [nn,X0,permuteDims,multiNetInfo] = ...
            aux_buildJointNetwork(modelPath,X0,vnnlibInfo,verbose);
        if verbose
            fprintf(' done\n');
        end
    else
        % ---- Single-network path -----------------------------------------
        if verbose
            fprintf('--- Loading network...');
        end
        [nn,permuteDims] = aux_loadNetwork(benchName_,modelPath,vnnlibPath,verbose);
        if verbose
            fprintf(' done\n');
        end
        multiNetInfo = [];
    end

    if verbose
        fprintf('--- Storing MATLAB file...');
    end
    % Create filename.
    instanceFilename = getInstanceFilename(benchName,modelPath,vnnlibPath);

    save(instanceFilename,'nn','options','permuteDims','X0','specs', ...
        'multiNetInfo','vnnlibInfo');
    if verbose
        fprintf(' done\n');
        aux_printOptions(options);
    end

catch e
    fprintf(newline);
    fprintf(e.message);
    fprintf(newline);
    res = 1;
    return;
end

res = 0;
end


% Auxiliary functions -----------------------------------------------------

function [nn,permuteDims] = aux_loadNetwork(benchName,modelPath,vnnlibPath,verbose)
permuteDims = false;
% model name needed for nn4sys instance filtering
[~,modelName,~] = getInstanceFilename(benchName,modelPath,vnnlibPath);

% Resolve .gz compressed ONNX (2026 benchmarks store compressed files).
% Multi-network paths (starting with '[') are handled by aux_buildJointNetwork (part B).
if ~startsWith(modelPath,'[')
    modelPath = aux_resolveVnncompPath(modelPath);
end

% VNN-COMP'24 Benchmarks ------------------------------------------------
if strcmp(benchName,'test') ...
        || strcmp(modelName{1},'test_nano') % is called after each benchmark
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'','');

elseif strcmp(benchName,'acasxu')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BSSC');

elseif strcmp(benchName,'cctsdb_yolo')
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'cgan')
    % --- TODO: implement convTranspose
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'cifar100')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS', ...
        '','dagnetwork',true);
    permuteDims = true;

elseif strcmp(benchName,'collins_aerospace_benchmark') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'collins_rul_cnn')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
    permuteDims = true;

elseif strcmp(benchName,'collins_yolo_robustness') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'cora')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'dist_shift')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'linearizenn')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC', ...
        '','dagnetwork',true);

elseif strcmp(benchName,'lsnc') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'metaroom')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
    permuteDims = true;

elseif strcmp(benchName,'ml4acopf') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'ml4acopf_2024') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'nn4sys')
    if ~strcmp(modelName{1},'lindex') && ...
            ~strcmp(modelName{1},'lindex_deep')
        throw(CORAerror('CORA:notSupported',...
            sprintf("Model '%s' of benchmark '%s' is not " + ...
            "supported!",modelPath,benchName)));
    end
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC','BC');

elseif strcmp(benchName,'safenlp')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'tinyimagenet')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS', ...
        '','dagnetwork',true);
    permuteDims = true;

elseif strcmp(benchName,'tllverifybench')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'traffic_signs_recognition') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'vggnet16') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'vit') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

elseif strcmp(benchName,'yolo') % not supported
    throw(CORAerror('CORA:notSupported',...
        sprintf("Benchmark '%s' not supported!",benchName)));

% VNN-COMP'22 Benchmarks ------------------------------------------------
elseif strcmp(benchName,'mnist_fc')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'SSC');

elseif strcmp(benchName,'oval21')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
    permuteDims = true;

elseif strcmp(benchName,'reach_prob_density')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');

elseif strcmp(benchName,'rl_benchmarks')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BC');
    % dubins rejoin instances not supported
    if ismember({vnnlibPath}, ...
            {'vnnlib/dubinsrejoin_case_safe_10.vnnlib', ...
            'vnnlib/dubinsrejoin_case_safe_13.vnnlib', ...
            'vnnlib/dubinsrejoin_case_safe_15.vnnlib', ...
            'vnnlib/dubinsrejoin_case_safe_16.vnnlib', ...
            'vnnlib/dubinsrejoin_case_safe_17.vnnlib'})
        throw(CORAerror('CORA:notSupported',...
            sprintf("Specification '%s' of benchmark '%s' is not " + ...
            "supported!",vnnlibPath,benchName)));
    end

elseif strcmp(benchName,'sri_resnet_a')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
    permuteDims = true;

elseif strcmp(benchName,'sri_resnet_b')
    nn = neuralNetwork.readONNXNetwork(modelPath,verbose,'BCSS');
    permuteDims = true;

else
    throw(CORAerror('CORA:notSupported',...
        sprintf("Unknown benchmark '%s'!",benchName)));
end

end

function gpu = aux_isGPUavailable()
try
    if ~isempty(which('gpuDeviceCount'))
        gpu = gpuDeviceCount('available') > 0;
    else
        gpu = false;
    end
catch
    gpu = false;
end
end

function [name,year] = aux_parseBenchmarkName(str) % strip optional year suffix
tok = regexp(str, '^(.+?)[_-](20\d{2})$', 'tokens');
if ~isempty(tok)
    name = tok{1}{1};
    year = tok{1}{2};
else
    name = str;
    year = '';
end
end

function [nn, X0_reduced, permuteDims, netInfo] = aux_buildJointNetwork(modelPath, X0_joint, info, verbose)
% parse the Python-list modelPath, load both ONNX sub-networks, and delegate
% to nnHelper.buildJointNetwork for the joint network + input projection

permuteDims = false;

% Parse "[('f', 'path1'), ('g', 'path2')]" into role/path pairs.
toks = regexp(modelPath, '\(''([^'']+)'',\s*''([^'']+)''\)', 'tokens');
if numel(toks) ~= 2
    throw(CORAerror('CORA:converterIssue', ...
        sprintf('Expected exactly 2 networks in modelPath, got %d', numel(toks))));
end
pathF = aux_resolveVnncompPath(toks{1}{2});
pathG = aux_resolveVnncompPath(toks{2}{2});

if verbose
    fprintf('\n  f: %s', pathF);
    fprintf('\n  g: %s\n', pathG);
end
nn_f = neuralNetwork.readONNXNetwork(pathF, verbose, 'BSSC');
nn_g = neuralNetwork.readONNXNetwork(pathG, verbose, 'BSSC');

% Assemble the joint network and project the input polytope (pure logic).
[nn, X0_reduced, netInfo] = nnHelper.buildJointNetwork(nn_f, nn_g, X0_joint, info);
end

function resolvedPath = aux_resolveVnncompPath(relPath)
% resolve a path that may be .gz-compressed: try as-is, then '.gz', then with
% one intermediate directory component stripped (e.g. onnx/original/ -> onnx/)

resolvedPath = relPath;
if exist(relPath,'file'), return; end

% Try with .gz suffix.
if exist([relPath '.gz'],'file')
    resolvedPath = aux_decompressGz([relPath '.gz']);
    return;
end

% Try stripping one intermediate path component (e.g. 'onnx/original/foo'
% -> 'onnx/foo').
parts = strsplit(relPath,'/');
if numel(parts) >= 3
    stripped = strjoin([parts(1), parts(3:end)],'/');
    if exist(stripped,'file')
        resolvedPath = stripped;
        return;
    end
    if exist([stripped '.gz'],'file')
        resolvedPath = aux_decompressGz([stripped '.gz']);
        return;
    end
end

% Return original path and let the caller produce the error.
end

function outPath = aux_decompressGz(gzPath)
% Decompress a .gz file next to the archive (idempotent).
outPath = gzPath(1:end-3);
if ~exist(outPath,'file')
    [outDir,~,~] = fileparts(outPath);
    if isempty(outDir), outDir = '.'; end
    gunzip(gzPath, outDir);
end
end

function aux_printOptions(options) % print options summary table
table = CORAtableParameters('neuralNetwork/verify options');
table.printHeader();
table.printContentRow('GPU',string(options.nn.train.use_gpu));
table.printContentRow('Poly. Method',options.nn.poly_method);
table.printContentRow('Batchsize', ...
    string(options.nn.train.mini_batch_size));
table.printContentRow('Interval Center', ...
    string(options.nn.interval_center));
% generator / approximation error settings
table.printContentRow('Num. init. Generators', ...
    string(options.nn.train.num_init_gens));
table.printContentRow('Num. approx. Error (per nonl. Layer)', ...
    string(options.nn.train.num_approx_err));
table.printContentRow('approx. Error Heuristic', ...
    options.nn.approx_error_order);
table.printContentRow('Falsification Method', ...
    options.nn.falsification_method);
table.printContentRow('Refinement Method', ...
    options.nn.refinement_method);
table.printContentRow('Num. of Splits', ...
    string(options.nn.num_splits));
table.printContentRow('Num. of Dimensions', ...
    string(options.nn.num_dimensions));
table.printContentRow('Num. of Neuron-Splits', ...
    string(options.nn.num_neuron_splits));
table.printContentRow('Num. of ReLU-Tightening Constraints', ...
    string(options.nn.num_relu_constraints));
table.printFooter();
end

% ------------------------------ END OF CODE ------------------------------
