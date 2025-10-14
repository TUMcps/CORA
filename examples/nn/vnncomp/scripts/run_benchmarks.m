function run_benchmarks(benchmarks,datapath,resultspath)
% run_benchmarks - run all benchmarks.
%
% Syntax:
%    run_benchmarks(benchmarks,datapath,resultspath)
%
% Inputs:
%    benchmarks - names of the benchmarks
%    datapath - path to the benchmark files
%    resultsPath - path to the results directory
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Run all benchmarks in the current directory.

% Restrict number of CPU threads.
maxNumCompThreads(4);

% Get the basepath.
basepath = pwd;
% Prepend the basepath to the results path.
resultspath = sprintf('%s/%s',basepath,resultspath);
% Create a results directory.
mkdir(resultspath);

% Specify base directory.
benchpaths = {...
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

for i=1:length(benchdirs)
    % Get the benchmark name.
    benchnamei = benchdirs(i).name;
    benchpathi = benchdirs(i).folder;
    fprintf('__________________________________________________________________\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('BENCHMARK %s (%d/%d)\n',benchnamei,i,length(benchdirs));
    fprintf('------------------------------------------------------------------\n');
    % Create a results directory.
    benchresultpath = sprintf('%s/2025_%s',resultspath,benchnamei); 
    mkdir(benchresultpath);
    % Change directory to the current benchmark.
    cd(sprintf('%s/%s',benchpathi,benchnamei));
    % Run all instances of the benchmark.
    [numVerif,numFals,numUnknown] = ...
        run_instances(benchnamei,benchresultpath);
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
