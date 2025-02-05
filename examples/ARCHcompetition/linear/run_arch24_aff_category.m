function run_arch24_aff_category()
% run_arch24_aff_category - runs all benchmarks from the linear category
%     from the 2024 ARCH competition
%
% Syntax:
%    text = run_arch24_aff_category()
%
% Inputs:
%    -
%
% Outputs:
%    -

% Authors:       Mark Wetzlinger
% Written:       23-March-2023
% Last update:   30-May-2024
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% regular directory
outputdir = [CORAROOT filesep 'examples' filesep 'ARCHcompetition' ...
    filesep 'linear' filesep 'results' filesep];
mkdir(outputdir)

% initialize .csv file
fid = fopen([outputdir 'results.csv'],'w');

% helper dashes
longdashes = '-------------------------------------------------------------------';
shortdashes = '------------------';

% header
text = 'benchmark,instance,result,time';
fprintf(fid,'%s\n',text);

% different benchmark groups (order matters!) ---
% name
benchmark_names = {'rand'; 'heat3D'; 'beam'; 'iss'; 'platoon'; 'rendezvous'; ...
    'powerTrain'; 'gearbox'; 'brake'};
% ids
benchmark_ids = { ...
    {'RAND01'; 'RAND02'}; ...
    {'HEAT01'; 'HEAT02'; 'HEAT03'}; ...
    {'CBC01'; 'CBC02'; 'CBC03'; 'CBF01'; 'CBF02'; 'CBF03'}; ...
    {'ISSC01_ISS02'; 'ISSC01_ISU02'; 'ISSF01_ISS01'; 'ISSF01_ISU01'}; ...
    {'PLAA01_BND42'; 'PLAA01_BND50'; 'PLAD01_BND30'; 'PLAD01_BND42'; 'PLAN01'}; ...
    {'SRA01'; 'SRA02'; 'SRA03'; 'SRA04'; 'SRA05'; 'SRA06'; 'SRA07'; 'SRA08'; ...
        'SRNA01'; 'SRU01'; 'SRU02'}; ...
    {'DTN01'; 'DTN02'; 'DTN03'; 'DTN04'; 'DTN05'; 'DTN06'}; ...
    {'GRBX01'; 'GRBX02'}; ...
    {'BRKDC01'; 'BRKNC01'; 'BRKNP01'}};
% respective functions
benchmark_funcs = {...
    'linear_verify'; ...
    {'linear_verifyFast'; 'linear_verifyFast'; 'linear_reach'}; ...
    'linear_verifyFast'; ...
    'linear_verifyFast'; ...
    'linearParam_reach'; ...
    'hybrid_reach'; ...
    'hybrid_reach'; ...
    'hybrid_reach'; ...
    'hybrid_reach'};
% keep plots
benchmark_plots = {
    [false; false]; ...
    [true; false; false]; ...
    [true; false; false; true; false; false]; ...
    [true; true; true; true]; ...
    [true; true; true; false; true]; ...
    [true; false; false; false; false; false; false; false; true; false; false]; ...
    [false; false; true; false; false; false]; ...
    [true; false]; ...
    [true; false; false]};
% build benchmark struct
benchmarks = struct(...
    'name', benchmark_names, ...
    'id', benchmark_ids, ...
    'prefix', benchmark_funcs, ...
    'plots', benchmark_plots);

% plots in ARCH23: HEAT01, CBC01, CBF01, ISSC, ISSF, SRNA01, SRA01, DTN03,
%                  PLAA01, PLAD01, PLAN01, GRBX01, BRKDC01

% counter for number of failed benchmark executions
cnt = 0;

% loop over all benchmark groups
for g=1:length(benchmarks)
    disp([benchmarks(g).name ' ' longdashes]);

    % loop over instance
    for i=1:length(benchmarks(g).id)
        disp([benchmarks(g).id{i} ' ' shortdashes]);
        
        % compose function name
        if iscell(benchmarks(g).prefix)
            algname = benchmarks(g).prefix{i};
        else
            algname = benchmarks(g).prefix;
        end
        % keep '23' in the name for now (renaming is time-consuming...)
        funcname = ['benchmark_' algname '_ARCH23_' benchmarks(g).name ...
            '_' benchmarks(g).id{i}];
        
        try 
            eval(['text = ' funcname '();']);
        catch ME
            disp(['Benchmark failed: ' ...
                benchmarks(g).name ', ' benchmarks(g).id{i}]);
            cnt = cnt + 1;
            continue
        end

        % print result to .csv file
        try
            fprintf(fid,'%s\n',text);
        catch ME
            disp(['Writing to .csv file failed: ' ...
                benchmarks(g).name ', ' benchmarks(g).id{i}]);
            continue
        end

        % save figure (if given)
        if benchmarks(g).plots(i)
            try
                saveas(gcf, [outputdir filesep benchmarks(g).id{i} '.png']);
                close;
            catch ME
                disp("Figure could not be saved/No figure to save.");
            end
        end
    end
    disp(" ");
end

% Close .csv file
try
    fclose(fid);
catch ME
    disp("File could not be closed!");
end

% final message
disp("Number of failed executions: " + cnt);
disp("--- Run over ---");

end

% ------------------------------ END OF CODE ------------------------------
