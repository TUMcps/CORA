function completed = main_vnncomp(varargin)
% main_vnncomp - runs all scripts to repeat the evaluation section of the paper
%
% results and plots will be saved to ./results
%
% Syntax:
%    completed = main_vnncomp()
%    completed = main_vnncomp(evalname)
%
% Inputs:
%    evalname - str, results are stored in ./results/<evalname>
%                    defaults to current date and time
%
% Outputs:
%    completed - boolean

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. SETTINGS (update as needed) ------------------------------------------

% general settings
PAPER_TITLE = 'CORA'; % !
VENUE_NAME = 'VNN-COMP';   % !
aux_runStartup(PAPER_TITLE,VENUE_NAME);

% also change 3. Scripts below

% 2. SETUP (nothing to change here) ---------------------------------------

% parse input
if nargin < 1
    evalname = 'cora';
else
    evalname = varargin{1};
end

% set up paths
basepath = '.';

% PATH                      VARIABLE        PURPOSE
% ./                        basepath        base path
% - ./code                  codepath        path to code
%   - ./cora                -               path to CORA
%   - ./scripts             -               path to auxiliary scripts
%   - ./main.m              -               main evaluation file (this file)
% - ./data                  datapath        path to data
% - ./results/<evalname>    resultspath     path to results
%   - ./evaluation          evalpath        path to store any evaluation results
%   - ./plots               plotspath       path to plots; any open figure after
%                                           each script will be stored there
%   - ./results.txt         -               logs all outputs to command window
%
[codepath,datapath,resultspath] = aux_setup(basepath,evalname);

% 3. RUN SCRIPTS (update as needed) ---------------------------------------

benchmarks = {...
    'test', ...
    ... % VNN-COMP'25 benchmarks
    ... 'acasxu_2023', ...
    ... 'cctsdb_yolo_2023', ... % (not supported; not main track)
    ... 'cersyve', ... % (test)
    ... 'cgan_2023', ... % (not supported)
    ... 'cifar100_2024', ... % (test & fix)
    ... 'collins_aerospace_benchmark', ... % (not supported)
    ... 'collins_rul_cnn_2022', ...
    ... 'cora_2024', ...
    ... 'dist_shift_2023', ...
    ... 'linearizenn_2024', ...
    ... 'lsnc_relu', ... % (not supported)
    ... 'malbeware', ...
    ... 'metaroom_2023', ...
    ... 'ml4acopf_2024', ... % (not supported)
    ... 'sat_relu', ...
    ... 'nn4sys', ... % (not supported; TODO: missing convolution 1D)
    ... 'relusplitter', ... % (test)
    ... 'safenlp_2024', ...
    ... 'soundnessbench', ... % (TODO: tune parameters)
    ... 'tinyimagenet_2024', ...
    ... 'tllverifybench_2023', ...
    ... 'traffic_signs_recognition_2023', ... % (not supported; not main track)
    ... 'vggnet16_2022', ... % (not supported; not main track)
    ... 'vit_2023', ... % (not supported)
    ... 'yolo_2023', ... % (not supported)
};

scripts = { ...
    @() run_benchmarks(benchmarks,datapath,resultspath), 'run_benchmarks';
};

% run scripts
aux_runScripts(scripts)

% 4. WRAP UP (nothing to change here) -------------------------------------

aux_wrapup()
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_runStartup(PAPER_TITLE,VENUE_NAME)
    rng(1)
    warning off

    % show startup block
    disp(' ')
    aux_seperateLine()
    disp(' ')
    disp('Repeatability Package')
    fprintf("Paper: %s\n", PAPER_TITLE)
    fprintf("Venue: %s\n", VENUE_NAME)
    fprintf("Date: %s\n", datestr(datetime()))
    disp(' ')
    if ~isempty(which('CORAVERSION'))
        fprintf('CORA: %s\n', CORAVERSION)
    end
    fprintf('Matlab: %s\n', version)
    fprintf('System: %s\n', computer)
    fprintf('GPU available: %i', canUseGPU)
    disp(' ')
    aux_seperateLine()
    disp(' ')
    pause(2) % to make the startup block readable
end

function [codepath,datapath,resultspath] = aux_setup(basepath,evalname)

    % set up paths
    codepath = sprintf("%s/code", basepath);
    datapath = sprintf("%s/data", basepath);
    resultspath = sprintf("%s/results/%s/%s", basepath, ...
        datestr(datetime,'yymmdd-hhMMss'), evalname);
    mkdir(resultspath)
    
    % for smooth images (only for eps)
    set(0, 'defaultFigureRenderer', 'painters')
    
    % set up diary
    resultstxt = sprintf("%s/results.txt", resultspath);
    delete(resultstxt)
    diary(resultstxt)
end

function aux_runScripts(scripts,plotspath,plotSettings)
    % run scripts

    % process input
    n = size(scripts, 1);
    fprintf("Running %d scripts.. \n", n);
    disp(" ")
    
    for i = 1:n
        % run script i
        aux_seperateLine()
        disp(' ')
        script = scripts{i, 1};
        name = scripts{i, 2};
    
        try
            % call script
            fprintf("Running '%s' ...\n", name)
            script();
            disp(" ")
            fprintf("'%s' was run successfully!\n", name)
    
        catch ME
            % error handling
            disp(" ")
            fprintf("An ERROR occured during execution of '%s':\n", name);
            disp(ME.getReport())
            disp("Continuing with next script..")
        end
    
        disp(" ")
    end
end

function aux_wrapup()
    % wrap up evaluation
    aux_seperateLine()
    disp(" ")
    disp("Completed!")
    fprintf("Date: %s\n", datestr(datetime()))
    diary off;
end

function aux_seperateLine()
    disp ------------------------------------------------------------------
end

% ------------------------------ END OF CODE ------------------------------
