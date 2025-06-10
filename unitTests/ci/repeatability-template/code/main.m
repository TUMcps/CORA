function completed = main(varargin)
% main - runs all scripts to repeat the evaluation section of the paper
%
% results and plots will be saved to ./results
%
% Syntax:
%    completed = main()
%    completed = main(evalname)
%
% Inputs:
%    evalname - str, results are stored in ./results/<evalname>
%                    defaults to current date and time
%
% Outputs:
%    completed - boolean

% ------------------------------ BEGIN CODE -------------------------------

% 1. SETTINGS (update as needed) ------------------------------------------

% general settings
PAPER_TITLE = '<paper-title>'; % !
VENUE_NAME = '<venue-name>';   % !
aux_runStartup(PAPER_TITLE,VENUE_NAME);

% 2. SETUP (nothing to change here) ---------------------------------------

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
[evalname,basepath,codepath,datapath,resultspath,evalpath,plotspath,plotSettings] = aux_setup(varargin{:});

% 3. RUN SCRIPTS (update as needed) ---------------------------------------

scripts = {; ...
    % list all evaluation scripts here:
    % "<scriptname>", @my_function_handle;
    % (<scriptname> is used for display and naming figures etc.)
    % e.g.
    "figure_1", @aux_mysimplefunc;
    "sec5_2", @() aux_mycomplexfunc(datapath,evalpath);
    };

% run scripts
aux_runScripts(scripts,plotspath,plotSettings)

% 4. WRAP UP (nothing to change here) -------------------------------------

aux_wrapup()
completed = true;

end


% SCRIPTS -----------------------------------------------------------------
% place all scripts here / call to external scripts in ./code/scripts

function aux_mysimplefunc()
    % re-create figure 1
    V = [ ...
     -2.000, 0.000, 2.000, 2.000, 0.000, -2.000 -2.000 ; ...
     -2.000, -2.000, 0.000, 2.000, 2.000, 0.000 -2.000 ; ...
     ];
    figure; hold on;
    plot(V(1,:),V(2,:));
end

function aux_mycomplexfunc(datapath, evalpath)

    % read stuff from data ---
    % ...

    % compute something
    % ...

    % obtain (dummy) results
    results = struct;
    results(1).val1 = 'v11';
    results(1).val2 = 'v12';
    results(2).val1 = 'v21';
    results(2).val2 = 'v22';
    results(3).val1 = 'v31';
    results(3).val2 = 'v32';

    % save raw results for later usage
    save(sprintf('%s/results.mat',evalpath),"results");

    % print table to command window 
    % (you might want to use CORAtable for that)
    for i=1:numel(results)
        disp(results(i))
    end

end


% Auxiliary functions -----------------------------------------------------

function aux_runStartup(PAPER_TITLE,VENUE_NAME)
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
    fprintf('GPU available: %i\n', canUseGPU)
    disp(' ')
    aux_seperateLine()
    disp(' ')
    pause(2) % to make the startup block readable
end

function [evalname,basepath,codepath,datapath,resultspath,evalpath,plotspath,plotSettings] = aux_setup(varargin)
    % determine evaluation name, set up paths, and define plot settings

    % parse input
    if nargin < 1
        evalname = datestr(datetime,'yymmdd-hhMMss');
    else
        evalname = varargin{1};
    end

    % set up paths
    basepath = '.';

    % set up paths
    codepath = sprintf("%s/code", basepath);
    datapath = sprintf("%s/data", basepath);
    resultspath = sprintf("%s/results/%s", basepath, evalname);
    mkdir(resultspath)
    evalpath = sprintf("%s/evaluation", resultspath);
    mkdir(evalpath)
    plotspath = sprintf("%s/plots", resultspath);
    mkdir(plotspath)
    
    % for smooth images (only for eps)
    set(0, 'defaultFigureRenderer', 'painters')
    
    % set up diary
    resultstxt = sprintf("%s/results.txt", resultspath);
    if exist("resultstxt","file")
        delete(resultstxt)
    end
    diary(resultstxt)

    % set up plotting settings
    plotSettings = struct;
    plotSettings.saveOpenFigures = true;
    plotSettings.saveAsFig = true;
    plotSettings.saveAsPng = true;
    plotSettings.saveAsEps = false;

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
        name = scripts{i, 1};
        script = scripts{i, 2};
    
        try
            % call script
            fprintf("Running '%s' ...\n", name)
            script();
            disp(" ")
            fprintf("'%s' was run successfully!\n", name)
    
            % save open figures
            if plotSettings.saveOpenFigures
                fprintf("Saving plots to '%s'..\n", plotspath)
                
                % find all open figures
                h = findobj('type', 'figure');
                m = length(h);
        
                % get unique name
                figure_name = sprintf('%s/%s', plotspath, name);

                % iterate through figures
                for j = 1:m
                    % save as desired extensions
                    if plotSettings.saveAsFig % .fig
                        savefig(sprintf("%s_%d.%s", figure_name, j, 'fig'));
                    end
                    if plotSettings.saveAsPng % .png
                        saveas(gcf, sprintf("%s_%d.%s", figure_name, j, 'png'));
                    end
                    if plotSettings.saveAsEps % .eps
                        saveas(gcf, sprintf("%s_%d.%s", figure_name, j, 'eps'), 'epsc');
                    end
                    % close figure
                    close(gcf)
                end
            end
    
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
