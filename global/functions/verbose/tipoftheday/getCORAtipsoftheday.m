function tips = getCORAtipsoftheday(rs)
% getCORAtipsoftheday - returns all CORA tips (of the day)
%
% Syntax:
%    tips = getCORAtipsoftheday
%    tips = getCORAtipsoftheday(rs)
%
% Inputs:
%    rs - random number generator
%
% Outputs:
%    tips - cell array with all tips
%

% Authors:       Tobias Ladner
% Written:       19-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 1
    % init random number generator
    rs = RandStream('mt19937ar', 'Seed', randi(2^32)); 
end

% 1 in 10 chance to check if new update is available
% (is quite slow due to url call, so not always computed)
tips = [];
if rs.rand() <= 0.1
    tips = aux_addNewRelease();
end

% gather static tips if previous step did not produce any tips
if isempty(tips)
    % gather tips by category
    tips = [
        % folder
        aux_getTipsApp()
        aux_getTipsContDynamics()
        aux_getTipsContSet()
        aux_getTipsConverter()
        aux_getTipsDiscDynamics()
        aux_getTipsExamples()
        aux_getTipsGlobal()
        aux_getTipsHybridDynamics()
        aux_getTipsManual()
        aux_getTipsMatrixSets()
        aux_getTipsModels()
        aux_getTipsNN()
        aux_getTipsSpecification()
        aux_getTipsUnitTests()
        % other
        aux_getTipsSolver()
        aux_getTipsPlot()
        % ___ of the day
        aux_getClassOfTheDayCandidates(rs, 10)
        aux_getContSetFunctionOfTheDayCandidates(rs, 10)
        aux_getExampleOfTheDayCandidates(rs, 10)
    ];
end

end


% Auxiliary functions -----------------------------------------------------

function tips = aux_getTipsApp()
    tips = {
        ['CORA provides a graphical user interface to verify dynamic systems.\n' ...
        'To start the app, type <a href="matlab:coraApp">coraApp</a> in the command window.\n' ...
        aux_addVisitManualText('Sec. 9')]
    };
end

function tips = aux_getTipsContDynamics()
    tips = {
        ['CORA computes reachable sets for linear systems, nonlinear systems as well as for systems with constraints. Continuous as well as discrete time models are supported.\n' ...
        'Uncertainty in the system inputs as well as uncertainty in the model parameters can be explicitly considered.\n' ...
        'In addition, CORA also provides capabilities for the simulation of dynamical models.' ...
        aux_addVisitManualText('Sec. 4')],
        ['Don''t know how to model your real-world system? \n' ...
        'Utilize our ' aux_addOpenFileText('contDynamics/conform') ' function to identify reachset-conformant models from data. \n' ...
        'Based on your prior system knowledge, you can choose to identify only the \n' ...
        'uncertainty sets (white-box identification) or extend your identification to \n' ...
        'unknown model parameters (gray-box identification) or the entire model structure \n' ...
        '(black-box identification).']
    };
end

function tips = aux_getTipsContSet()
    tips = {
        ['CORA has a modular design, making it possible to use the capabilities of the various set representations for other purposes besides reachability analysis.\n' ...
        'The toolbox implements vector set representation, e.g., intervals, zonotopes, Taylor models, and polytopes, \n' ...
        'as well as matrix set representations such as matrix zonotopes and interval matrices.\n' ...
        aux_addVisitManualText('Sec. 2')],
        ['Computed sets in the workspace can be printed into the command window for easy re-initialization.\n' ...
            aux_addHelpCWText('printSet')]
    };
end

function tips = aux_getTipsConverter()
    tips = {
        ['CORA provides capabilities to convert from and to various file formats.\n' ...
        'This includes SpaceEx, CommonRoad, CommonOcean, and ONNX Neural Networks.\n' ...
        aux_addVisitManualText('Sec. 8')]
    };
end

function tips = aux_getTipsDiscDynamics()
    tips = {
        ['In addition to continuous-time systems, CORA also supports discrete-time systems.\n' ...
        aux_addVisitManualText('Sec. 4.2')]
    };
end

function tips = aux_getTipsExamples()
    tips = {
        ['CORA provides a variety of examples to give you hands-on experience in CORA.\n' ...
        'Visit the ./examples folder for more information.']
    };
end

function tips = aux_getTipsGlobal()
    tips = {
        ['Do you need to update Matlab? Use our installation script to quickly get you back to speed.\n' ...
            aux_addHelpCWText('installCORA')]
        ['Do you require CORA in a repeatability package? We got you covered using docker.\n' ...
        sprintf('Check out %s and the accompanying %s', ...
        aux_addOpenFileText('./cora/unitTests/ci/build/Dockerfile'), ...
        aux_addOpenFileText('./cora/unitTests/ci/build/README.md')) '.']
        ['CORA can accumulate converted system dynamics and auxiliary files over time. Consider resetting CORA regularly.\n' ...
            aux_addHelpCWText('resetCORA')]
        ['Do you want to share a reachable set you obtained with CORA?\n' ...
        'Consider creating a video of the reachable set evolving over time using the ' aux_addOpenFileText('recordCORAvideo') ' function.\n' ...
        'An example can be found here: ' aux_addLinkText('https://www.youtube.com/watch?v=BqRy4bBGn6E')]
        ['Today we want to guide you to our website: ' aux_addLinkText('https://cora.in.tum.de') ',\n' ...
        'where you can find examples, older CORA versions, and publication using CORA.\n' ...
        'Did you use CORA in one of your publications? Tell us and we are happy to add it to our website.']
        ['It is advisable to consider numerical issues when writing code.\n' ...
        'CORA provides a tolerance option for most functions to compensate for small numerical issues.\n' ...
        aux_addHelpCWText({'withinTol','compareMatrices'})]
        ['CORA is continuously developed with bug fixes and new features.\n' ...
        'To make sure that all folders of CORA are on the Matlab path,\n' ...
        'we encourage you to add ' aux_addOpenFileText('updateCORApath') ' to your ' aux_addOpenFileText('startup') ' file.']
        ['Do you want to display some information in a nice format?\n' ...
        aux_addOpenFileText('CORAtable') 's offers you multiple design options.\n' ...
        aux_addHelpCWText('CORAtable')]

        % macros
        ['CORA performs rigorous input argument checking to ensure correctness and provide valuable error message to the user.\n' ...
        'However, these can be quite slow once your code is up and running.\n' ...
        'You can disable the checks via the ' aux_addOpenFileText('CHECKS_ENABLED')  ' macro.']
        ['Do you need to dynamically find the root folder of CORA, e.g., for repeatability?\n' ...
        'The ' aux_addOpenFileText('CORAROOT') ' macro does just that.\n' ...
        aux_addHelpCWText('CORAROOT')]

        % verification
        ['CORA provides a variety of STL verification methods based on reachable sets.\n' ...
        'Have a look at ' aux_addOpenFileText('modelChecking') ' to see the available algorithms.']
        ['You can use CORA to verify STL properties for hybrid systems.\n' ...
        'See ' aux_addOpenFileText('example_stl_bouncingBall') ' for an example.']
    };
end

function tips = aux_getTipsHybridDynamics()
    tips = {
        ['CORA is capable to calculate the reachable sets for hybrid systems.\n' ...
        'All implemented dynamic system classes can be used to describe the different continuous flows for the discrete system states.\n' ...
        'Further, multiple different methods for the calculation of the intersections with guard sets are implemented in CORA.\n' ...
        aux_addVisitManualText('Sec. 4.3')]
    };
end

function tips = aux_getTipsManual()
    tips = {
        ['We document the capabilities of CORA in the <a href="cora.in.tum.de/manual">CORA manual</a>.\n' ...
        'Please check the manual for definitions, details, examples, and more.']
    };
end

function tips = aux_getTipsMatrixSets()
    tips = {
        ['CORA also has some matrix set representations implemented.\n' ...
        aux_addVisitManualText('Sec. 3.2')]
    };
end

function tips = aux_getTipsModels()
    tips = {
        'All models used in CORA can be found in the respective folder in ./models.'
    };
end

function tips = aux_getTipsNN()
    tips = {
        ['CORA enables the formal verification of neural networks, both in open-loop as well as in closed loop scenarios.\n' ...
        'Open-loop verification refers to the task where properties of the output set of a neural network are verified, e.g. correctly classified images given noisy input.\n' ...
        'In closed-loop scenarios, the neural network is used as controller of a dynamic system and is neatly integrated in the reachability algorithms above, e.g. controlling a car while keeping a safe distance.\n' ...
        aux_addVisitManualText('Sec. 6')]
        ['CORA enables the training of verifiably robust neural networks in the supervised and reinforcement learning setting.\n' ...
        'To get started, check out the example files ' aux_addOpenFileText('example_neuralNetwork_train') ' and ' aux_addOpenFileText('example_neuralNetwork_rl_agentDDPG_Quad1D') '.\n' ...
        aux_addVisitManualText('Sec. 6.4')]
    };
end

function tips = aux_getTipsSpecification()
    tips = {
        ['To define the specifications of your system, please check the ./specification folder.\n' ...
        'CORA enables users the specifications of an unsafe and safe sets using the @specification class, as well as STL and RTL.']
    };
end

function tips = aux_getTipsUnitTests()
    tips = {
        ['CORA is actively developed and continuously adds new functionalities as well as bug fixes.\n' ...
        'Make sure to regularily check if there is a new version.\n' ...
        'Found a bug? Please report it on <a href="https://github.com/TUMcps/CORA">GitHub</a>.']    
    };
end

function tips = aux_getTipsSolver()
    tips = {
        ['CORA is using third-party solvers for the computation of some functionalities.\n' ...
        'We primarily use YALMIP for solving optimization problems but you might want to install MOSEK to speed up certain computations.\n' ...
        aux_addVisitManualText('Sec. 1.3')]
    };
end

function tips = aux_getTipsPlot()
    tips = {
        ['CORA comes with various plotting capabilities of continuous sets.\n' ...
        aux_addHelpCWText('contSet/plot')]
        ['CORA enables the plotting of the computed reachable sets both over two dimenions as well as over time.\n' ...
        aux_addHelpCWText({'reachSet/plot', 'reachSet/plotOverTime'})]
        ['Plotting large reachable sets can take time.\n' ...
        'To speed up the process, consider unifying the sets with the ''Unify'' option.\n' ...
        aux_addHelpCWText({'reachSet/plot', 'reachSet/plotOverTime'})]
        ['Tired of using the same colors over and over?\n' ...
        'CORA enables customizing your plots do your needs.\n' ...
        sprintf('Visit %s for more information.', aux_addOpenFileText('example_plot_color'))]
        ['Do you want the CORA look-and-feel for your plots?\n' ...
        'Consider calling "useCORAcolors(''CORA:contDynamics'')" before plotting your reachable sets.\n' ...
        'This function also enables you to plot multiple reachable sets with distinct reachSet colors for comparison.\n' ...
        sprintf('Visit %s for an example.', aux_addOpenFileText('example_linear_reach_04_adaptive'))]
    };
end

function tips = aux_addNewRelease()
    tips = {};
    try
        % get latest release
        tags = webread("https://api.github.com/repos/TUMcps/CORA/tags");
        
        % check if github ping was successful
        if ~isempty(tags) && isstruct(tags) && isfield(tags(1),'name')
            % read version on github ('v202x.x.x')
            GITHUBVERSION = tags(1).name;
            % check if matches CORAVERSION ('CORA v202x.x.x')
            if ~contains(CORAVERSION,GITHUBVERSION)
                % add tip about new CORAVERSION
                tips{end+1} = sprintf( ...
                    ['A new CORA version is available: <strong>CORA %s</strong> (currently %s is installed).\\n' ...
                    'Make sure to update it soon to get the latest features and bug fixes!\\n' ...
                    'Read more about it here: ' aux_addLinkText('https://cora.in.tum.de/pages/archive#release-notes')], ...
                    GITHUBVERSION, CORAVERSION);
            end
        end
    catch ME
        % catch exception; output warning
        CORAwarning('CORA:global','Unable to check if a new CORA version is available.')
    end
end

% helper texts ---

function visitmanualtext = aux_addVisitManualText(sec)
    if nargin == 0
        visitmanualtext = 'Please visit the <a href="https://cora.in.tum.de/manual">CORA manual</a> for more information.';
    else
        visitmanualtext = sprintf('Please visit %s in the <a href="https://cora.in.tum.de/manual">CORA manual</a> for more information.',sec);
    end
end

function helpcwtext = aux_addHelpCWText(name)
    % convert to cell (for multiple texts)
    if ~iscell(name)
        name = {name};
    end

    % add link
    links = cellfun(@(name) sprintf('<a href="matlab:help %s">help %s</a>', name, name), name, 'UniformOutput',false);

    % join text
    linktext = strjoin(links, ' and ');

    helpcwtext = sprintf('For more information, type %s in the command window.', linktext);
end

function helpcwtext = aux_addDocText(name)
    helpcwtext = sprintf('For more information, type <a href="matlab:doc %s">doc %s</a> in the command window.', name, name);
end

function openfiletext = aux_addOpenFileText(file)
    openfiletext = sprintf('<a href="matlab:open %s">%s</a>', file, file);
end

function linktext = aux_addLinkText(link, linktext)
    if nargin < 2
        linktext = link;
    end
    linktext = sprintf('<a href="%s">%s</a>', link, linktext);
end

% '___ of the day' candidates ---

function tips = aux_getClassOfTheDayCandidates(rs, n)
    
    % find classes of cora
    coraroot = CORAROOT;
    classes = [
        dir([coraroot '/contDynamics/@*']);
        dir([coraroot '/contSet/@*']);
        % dir([coraroot '/discrDynamics/@*']);
        dir([coraroot '/hybridDynamics/@*']);
        dir([coraroot '/global/classes/@*']);
        dir([coraroot '/nn/@*']);
        dir([coraroot '/specification/@*']);
    ];

    % select n
    classes = classes(rs.randperm(numel(classes),n));

    % build tips
    tips = arrayfun(@(class) [ ...
        'Class of the day: ' ...
        aux_addOpenFileText(sprintf('%s', class.name(2:end))) '.\n' ...
        aux_addDocText(class.name(2:end))], ...
        classes, 'UniformOutput', false);
end

function tips = aux_getContSetFunctionOfTheDayCandidates(rs, n)
    
    % find functions of contSet
    funcs = dir([CORAROOT filesep 'contSet/@contSet/*.m']);

    % filter certain functions
    idx = arrayfun(@(file) ...
        ~endsWith(file.name,'_.m') ...
        && ~endsWith(file.name,'contSet.m') ...
        && ~endsWith(file.name,'display.m'), ...
        funcs, 'UniformOutput',true);
    funcs = funcs(idx);

    % select n
    funcs = funcs(rs.randperm(numel(funcs),n));

    % build tips
    tips = arrayfun(@(func) [ ...
        'Function of the day: ' ...
        aux_addOpenFileText(sprintf('contSet/%s', func.name(1:end-2))) '.\n' ...
        aux_addHelpCWText(['contSet/' func.name(1:end-2)])], ...
        funcs, 'UniformOutput', false);
end

function tips = aux_getExampleOfTheDayCandidates(rs, n)
    
    % find classes of cora
    examples = dir([CORAROOT filesep 'examples/**/*.m']);

    % select n
    examples = examples(rs.randperm(numel(examples),n));

    % build tips
    tips = arrayfun(@(example) [ ...
        sprintf('Example of the day: %s.', aux_addOpenFileText(example.name(1:end-2))) '\n' ...
        'Checking out existing examples is a good way to get hands-on experience in CORA.'], ...
        examples, 'UniformOutput', false);
end

% ------------------------------ END OF CODE ------------------------------
