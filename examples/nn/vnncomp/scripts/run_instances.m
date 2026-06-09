function [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath,options,instanceIds,timeoutMultiplier)
% run_instances - run all instances of a benchmark.
%
% Syntax:
%    [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath)
%    [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath,options)
%    [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath,options,instanceIds)
%    [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath,options,instanceIds,timeoutMultiplier)
%
% Inputs:
%    benchname - name of the benchmark
%    resultsPath - path to the results directory
%    options - (optional) fully resolved options struct. When omitted,
%              getDefaultVNNCOMPoptions(benchname) is used.
%    instanceIds - (optional) array of instance indices to run (default: all)
%    timeoutMultiplier - (optional) fraction of official timeout to use (default: 1.0)
%
% Outputs:
%    numVerif - number of verified instances
%    numFals - number of falsified instances
%    numUnknown - number of unknown instances
%
% References:
%    [1] VNN-COMP'24
%    [2] VNN-COMP'25
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   06-June-2026 (BK, multi-network path support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% When no options provided, load benchmark defaults.
if nargin < 3 || isempty(options)
    options = getDefaultVNNCOMPoptions(benchname);
end
if nargin < 5 || isempty(timeoutMultiplier)
    timeoutMultiplier = 1.0;
end

% Run all instances found in the current directory.
verbose = true;

% Obtain all instances (use fullfile(pwd,...) so MATLAB path is not searched).
filename = fullfile(pwd,'instances.csv');
instances = readtable(filename,'Delimiter',',');
% Rename columns
instances.Properties.VariableNames = {'model','vnnlib','timeout'};
% Obtain number of instances.
N = size(instances,1);

% Init results.
benchnames = {};
models = {};
vnnlibs = {};
prepTimes = [];
results = {};
verifTimes = [];
totalTimes = [];
numSubproblems = [];

% Count number of verified, falsified, and unknown instances.
numVerif = 0;
numFals = 0;
numUnknown = 0;

% Count number of patches needed to verify an instance.
numVerifPatches = 0;

% Use all instances if not specified.
if nargin < 4 || isempty(instanceIds)
    instanceIds = 1:N;
end

for i=instanceIds
    fprintf('__________________________________________________________________\n');
    fprintf('INSTANCE (%d/%d)\n',i,N);
    fprintf('------------------------------------------------------------------\n');
    % Extract current instance.
    instance = instances(i,:);
    modelPath = instance.model{1};
    vnnlibPath = instance.vnnlib{1};
    timeout = instance.timeout * timeoutMultiplier;

    % Create instance filename (match .onnx stem, excluding Python-list chars
    % so multi-network paths also work).
    modelName = regexp(modelPath,'([^/''"() ]+)(?=\.onnx)','match');

    if strcmp(benchname,'safenlp')
        if contains(modelPath,'medical')
            modelName{1} = ['medical_' modelName{1}];
        elseif contains(modelPath,'ruarobot')
            modelName{1} = ['ruarobot_' modelName{1}];
        end
    end

    vnnlibName = regexp(vnnlibPath,'([^/]+)(?=\.vnnlib)','match');
    instanceFilename = sprintf('%s/%s_%s.counterexample',...
        resultsPath,modelName{1},vnnlibName{1});

    % Prepare the current instance with fully resolved options.
    prepare_instance(benchname,modelPath,vnnlibPath,verbose,options);

    totalTime = tic;

    % Run the current instance.
    [resStr,res] = run_instance(benchname,modelPath,vnnlibPath,...
        instanceFilename,timeout,verbose);

    instanceTime = toc(totalTime);

    if strcmp(resStr,'unsat') || strcmp(resStr,'unknown')
        % There is no counterexample; delete the file.
        delete(instanceFilename);
    end

    % Store outputs.
    benchnames = [benchnames; benchname];
    models = [models; fullfile(benchname, modelPath)];
    vnnlibs = [vnnlibs; fullfile(benchname, vnnlibPath)];
    prepTimes = [prepTimes; res.prepTime];
    results = [results; resStr];
    % Use verify-only time from verify.m; fall back to instanceTime on error.
    if res.time >= 0
        verifTimes = [verifTimes; res.time];
    else
        verifTimes = [verifTimes; instanceTime];
    end
    totalTimes = [totalTimes; res.totalTime];
    numSubproblems = [numSubproblems; res.numSubproblems];

    % Increment counters.
    numVerif = numVerif + strcmp(resStr,'unsat');
    numFals = numFals + strcmp(resStr,'sat');
    numUnknown = numUnknown + strcmp(resStr,'unknown');

    numVerifPatches = numVerifPatches + res.numSubproblems;

    fprintf('------------------------------------------------------------------\n');
    fprintf('__________________________________________________________________\n');
end

% Print stats.
statsTable = CORAtableParameters('Stats');
statsTable.printHeader();
statsTable.printContentRow('avg. #Verified Branches',...
    string(numVerifPatches/length(instanceIds)));
statsTable.printContentRow('avg. Time',...
    string(sum(verifTimes)/length(instanceIds)));
statsTable.printFooter();

% Generate results table.
resultsTable = table(benchnames,models,vnnlibs,prepTimes,results, ...
    verifTimes,totalTimes,numSubproblems);
% Write to file.
writetable(resultsTable,sprintf('%s/results.csv',resultsPath));

end

% ------------------------------ END OF CODE ------------------------------
