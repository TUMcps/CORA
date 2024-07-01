function run_arch24_nln_category
% run_arch24_nln_category - runs all benchmarks from the nonlinear category
%     from the 2024 ARCH competition
%
% Syntax:
%    text = run_arch24_nln_category()
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

% local
outputdir = [CORAROOT filesep 'examples' filesep 'ARCHcompetition' ...
    filesep 'nonlinear' filesep 'results' filesep];
fid = fopen([outputdir 'results.csv'],'w');

% header
text = 'benchmark,instance,result,time,accuracy,timesteps';
fprintf(fid,'%s\n',text);

longdashes = '-------------------------------------------------------------------';

benchmark_names = {'Robertson'; 'vanDerPol'; 'laubloomis'; ...
    'lotkaVolterra'; 'spacecraft'; 'trafficscenario'};
benchmark_funcs = {'nonlinear'; 'nonlinear'; 'nonlinear'; ...
    'hybrid'; 'hybrid'; 'nonlinear'};
benchmarks = struct(...
    'name', benchmark_names, ...
    'prefix', benchmark_funcs);

cnt = 0;
% Robertson has an issue... abstractionError_adaptive
% all others are ok

% loop over benchmarks
for g=1:length(benchmarks)
    disp([benchmarks(g).name ' ' longdashes]);

    % compose name
    % keep '23' in the name for now (renaming is time-consuming...)
    funcname = ['benchmark_' benchmarks(g).prefix '_reach_ARCH23_' ...
        benchmarks(g).name];
    
    try 
        eval(['text = ' funcname '();']);
    catch ME
        disp(['Benchmark failed: ' ...
            benchmarks(g).name ', ' benchmarks(g).name]);
        cnt = cnt + 1;
    end

    % print result to .csv file
    try
        if iscell(text)
            for i=1:length(text)
                fprintf(fid,'%s\n',text{i});
            end
        else
            fprintf(fid,'%s\n',text);
        end
    catch ME
        disp(['Writing to .csv file failed: ' ...
            benchmarks(g).name ', ' benchmarks(g).name]);
    end

    % save figure (if given)
    try
        saveas(gcf, [outputdir filesep benchmarks(g).name '.png']);
        close;
    catch ME
        disp("Figure could not be saved/No figure to save.");
    end

    disp(' ');
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
