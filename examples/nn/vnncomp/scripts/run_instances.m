function [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath)
% run_instances - run all instances of a benchmark.
%
% Syntax:
%    [numVerif,numFals,numUnknown] = run_instances(benchname,resultsPath)
%
% Inputs:
%    benchname - name of the benchmark
%    resultsPath - path to the results directory
%
% Outputs:
%    numVerif - number of verified instances
%    numFals - number of falsified instances
%    numUnknown - number of unknown instances
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Run all instances found in the current directory.
verbose = true;

% Obtain all instances.
filename = 'instances.csv';
instances = readtable(filename,'Delimiter',',');
% Rename columns
instances.Properties.VariableNames = {'model','vnnlib','timeout'};
% Obtain number of instances.
N = size(instances,1);

% Init results.
benchnames = {};
models = {};
vnnlibs = {};
prepTimes = {};
results = {};
verifTimes = {};

% Count number of verified, falsified, and unknown instances.
numVerif = 0;
numFals = 0;
numUnknown = 0;

% Count number of patches needed to verify an instance.
numVerifPatches = 0;

% Specify the instance ids.
instanceIds = 1:N; % All
% instanceIds = 1:45; % acaxXu prop1
% instanceIds = 46:90; % acaxXu prop2
% instanceIds = 91:135; % acaxXu prop3
% instanceIds = 136:180; % acaxXu prop4

for i=instanceIds
    fprintf('__________________________________________________________________\n');
    fprintf('INSTANCE (%d/%d)\n',i,N);
    fprintf('------------------------------------------------------------------\n');
    % Extract current instance.
    instance = instances(i,:);
    modelPath = instance.model{1};
    vnnlibPath = instance.vnnlib{1};
    timeout = instance.timeout;

    % Create instance filename.
    modelName = regexp(modelPath,'([^/]+)(?=\.onnx$)','match');

    if strcmp(benchname,'safenlp')
        if contains(modelPath,'medical')
            modelName{1} = ['medical_' modelName{1}];
        elseif contains(modelPath,'ruarobot')
            modelName{1} = ['ruarobot_' modelName{1}];
        end
    end

    vnnlibName = regexp(vnnlibPath,'([^/]+)(?=\.vnnlib$)','match');
    instanceFilename = sprintf('%s/%s_%s.counterexample',...
        resultsPath,modelName{1},vnnlibName{1});

    % Prepare the current instance.
    prepare_instance(benchname,modelPath,vnnlibPath);

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
    models = [models; 'vnncomp2024_benchmarks/benchmarks/' ...
        benchname '/' modelPath];
    vnnlibs = [vnnlibs; 'vnncomp2024_benchmarks/benchmarks/' ...
        benchname '/' vnnlibPath];
    prepTimes = [prepTimes; 0];
    results = [results; resStr];
    verifTimes = [verifTimes; instanceTime];

    % Increment counters.
    numVerif = numVerif + strcmp(resStr,'unsat');
    numFals = numFals + strcmp(resStr,'sat');
    numUnknown = numUnknown + strcmp(resStr,'unknown');

    if isfield(res,'numVerified')
        numVerifPatches = numVerifPatches + res.numVerified;
    end

    fprintf('------------------------------------------------------------------\n');
    fprintf('__________________________________________________________________\n');
end

% Print stats.
statsTable = CORAtableParameters('Stats');
statsTable.printHeader();
statsTable.printContentRow('avg. #Verified Branches',...
    string(numVerifPatches/length(instanceIds)));
statsTable.printContentRow('avg. Time',...
    string(sum([verifTimes{:}])/length(instanceIds)));
statsTable.printFooter();

% Generate results table.
resultsTable = table(benchnames,models,vnnlibs,prepTimes,results, ...
    verifTimes);
% Write to file.
writetable(resultsTable,sprintf('%s/results.csv',resultsPath),...
    'WriteVariableNames',false);

end

% ------------------------------ END OF CODE ------------------------------
