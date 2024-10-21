function gp = gpinit(gp)
%GPINIT Initialises a run.
%
%   GP = GPINIT(GP)
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPCHECK, GPDEFAULTS, GPINITPARALLEL

%determine status of symbolic, parallel and stats toolboxes
[gp.info.toolbox.symbolic, gp.info.toolbox.parallel, gp.info.toolbox.stats] = gptoolboxcheck;

%process function nodes before run
gp = procfuncnodes(gp);

%throw an error if there are no inputs, p_ERC=0 and there are no arity zero
%functions active
if  gp.nodes.inputs.num_inp == 0 && gp.nodes.const.p_ERC == 0 && ...
        isempty(find(gp.nodes.functions.arity(logical(gp.nodes.functions.active)) == 0, 1))
    error('No terminals (inputs, constants or zero arity functions) have been defined for this run.');
end

%initialise some state and tracker variables
gp.state.count = 1;
gp.state.best.fitness = [];
gp.state.best.individual = [];
gp.state.run_completed = false;
gp.state.current_individual = [];
gp.state.std_devfitness = [];
gp.state.terminate = false;
gp.state.force_compute_theta = false;
gp.fitness.returnvalues = cell(gp.runcontrol.pop_size,1);

%process mutation probabilities vector
gp.operators.mutation.cumsum_mutate_par = cumsum(gp.operators.mutation.mutate_par);

%init. history variables
gp.results.history.bestfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.meanfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.std_devfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.about = 'Fitness on training data';
gp.results.history = orderfields(gp.results.history);

%best of run (on training data) fields
gp.results.best.fitness = [];
gp.results.best.individual = [];
gp.results.best.returnvalues = [];
gp.results.best.foundatgen = [];
gp.results.best.about = 'Best individual on training data';
gp.results.best = orderfields(gp.results.best);

%assign field holding fitnesses to gp structure
gp.fitness.values = zeros(gp.runcontrol.pop_size,1);
gp.fitness.complexity = zeros(gp.runcontrol.pop_size,1);

if strncmpi(func2str(gp.fitness.fitfun),'regressmulti',12)  && ~isfield(gp.userdata,'bootSampleSize')
    gp.userdata.bootSampleSize = size(gp.userdata.ytrain,1);   
end

%cache init
if gp.runcontrol.usecache
    gp = initcache(gp);
end

if ~gp.runcontrol.quiet
    fns = [];
    for i=1:length(gp.nodes.functions.active_name_UC);
        fns = [fns ' ' gp.nodes.functions.active_name_UC{i}];
    end
    
    if gp.selection.tournament.lex_pressure
        lex_inf = 'True';
    else
        lex_inf = 'False';
    end
    
    disp(' ');
    % disp('-------------------------------------------------------------------------');
    % disp('GPTIPS 2');
    % disp('Symbolic data mining platform for MATLAB');
    % disp('Copyright (C) Dominic Searson 2009-2015');
    % disp(' ');
    % disp('Contact: searson@gmail.com');
    % disp(' ');
    % 
    % disp('This program is free software: you can redistribute it and/or modify');
    % disp('it under the terms of the GNU General Public License as published by');
    % disp('the Free Software Foundation, either version 3 of the License, or');
    % disp('(at your option) any later version.');
    % disp(' ');
    % disp('This program is distributed in the hope that it will be useful,');
    % disp('but WITHOUT ANY WARRANTY; without even the implied warranty of');
    % disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
    % disp('GNU General Public License for more details: http://www.gnu.org/licenses');
    % disp(' ');
    disp('-------------------------------------------------------------------------- ');
    disp(' ');
    disp('Run parameters');
    disp('--------------');
    disp(['Population size:         ' int2str(gp.runcontrol.pop_size)]);
    disp(['Number of generations:   ' int2str(gp.runcontrol.num_gen)]);
    disp(['Number of runs:          ' int2str(gp.runcontrol.runs)]);
    
    if gp.runcontrol.parallel.auto
        disp('Parallel mode :          auto ');
    elseif gp.runcontrol.parallel.enable
        disp('Parallel mode :          manual ');
    else
        disp('Parallel mode :          off ');
    end
    
    if gp.selection.tournament.p_pareto == 0
        disp('Tournament type:         regular');
    else
        disp(['Tournament type:         Pareto (probability = ' ...
            num2str(gp.selection.tournament.p_pareto) ')']);
    end
    
    disp(['Tournament size:         ' int2str(gp.selection.tournament.size)]);
    disp(['Elite fraction:          ' num2str(gp.selection.elite_fraction)]);
    
    if gp.runcontrol.usecache
        disp('Fitness cache:           enabled');
    else
        disp('Fitness cache:           disabled');
    end
    
    disp(['Lexicographic selection: ' lex_inf]);
    disp(['Max tree depth:          ' int2str(gp.treedef.max_depth)]);
    disp(['Max nodes per tree:      ' int2str(gp.treedef.max_nodes)]);
    disp(['Using function set:     ' fns]);
    disp(['Number of inputs:        ' int2str(gp.nodes.inputs.num_inp)]);
    
    if gp.genes.multigene
        disp(['Max genes:               ' int2str(gp.genes.max_genes)]);
    end
    
    if ~gp.nodes.const.p_ERC
        disp('Using no constants');
    else
        disp(['Constants range:         [' num2str(gp.nodes.const.range) ']']);
    end
    
    if gp.fitness.complexityMeasure
        disp('Complexity measure:      expressional');
    else
        disp('Complexity measure:      node count');
    end
    
    disp(['Fitness function:        ' func2str(gp.fitness.fitfun) '.m']);
    disp(' ');
end

%parallel computing initialisation
gp = gpinitparallel(gp);

%log run start time
gp.info.startTime = datestr(now,0);

%set elapsed run time count to zero seconds
gp.state.runTimeElapsed = 0;

function gp=procfuncnodes(gp)
%PROCFUNCNODES Process required function node information prior to a run.

%loop through function nodes and generate arity list
for i = 1:length(gp.nodes.functions.name)
    
    arity = nargin(gp.nodes.functions.name{i});
    
    %some functions have a variable number of input arguments (e.g. rand)
    %In this case generate an error message and exit
    if arity == -1
        error(['The function ' gp.nodes.functions.name{i} ...
            ' may not be used (directly) as a function node because it has a variable number of arguments.']);
    end
    
    gp.nodes.functions.arity(i) = arity;
end

if ~isfield(gp.nodes.functions,'active') || isempty(gp.nodes.functions.active)
    gp.nodes.functions.active = ones(1,length(gp.nodes.functions.name));
end

gp.nodes.functions.active = logical(gp.nodes.functions.active);

%check max number of allowed functions not exceeded
gp.nodes.functions.num_active = numel(find(gp.nodes.functions.active));
if gp.nodes.functions.num_active > 22
    error('Maximum number of active functions allowed is 22');
end

%Generate single char Active Function IDentifiers (afid)(a->z excluding
%x,e,i,j) to stand in for function names whilst processing expressions.
%Exclusions are because 'x' is reserved for input nodes, 'e' is used for
%expressing numbers in standard form by Matlab and, by default, 'i' and 'j'
%represent sqrt(-1).
charnum = 96; skip = 0;
for i=1:gp.nodes.functions.num_active
    while true      %e                          %i                    %j                         %x
        if (charnum+i+skip)==101 || (charnum+i+skip)==105 || (charnum+i+skip)==106 || (charnum+i+skip)==120
            skip = skip + 1;
        else
            break
        end
        
    end
    afid(i) = char(charnum+i+skip);
end

%extract upper case active function names for later use
gp.nodes.functions.afid = afid;
temp = cell(gp.nodes.functions.num_active,1);

if numel(gp.nodes.functions.name) ~= numel(gp.nodes.functions.active)
    error('There must be the same number of entries in gp.nodes.functions.name and gp.nodes.functions.active. Check your config file.');
end

[temp{:}] = deal(gp.nodes.functions.name{gp.nodes.functions.active});
[gp.nodes.functions.active_name_UC] = upper(temp);

%generate index locators for arity >0 and arity == 0 active functions. The
%treegen function needs his info later for identifying which functions are
%terminal and which are internal nodes.
active_ar = (gp.nodes.functions.arity(gp.nodes.functions.active));
fun_argt0 = active_ar > 0;
fun_areq0 =~ fun_argt0;

gp.nodes.functions.afid_argt0 = gp.nodes.functions.afid(fun_argt0); %functions with arity > 0
gp.nodes.functions.afid_areq0 = gp.nodes.functions.afid(fun_areq0); %functions with arity == 0
gp.nodes.functions.arity_argt0 = active_ar(fun_argt0);

gp.nodes.functions.fun_lengthargt0 = numel(gp.nodes.functions.afid_argt0);
gp.nodes.functions.fun_lengthareq0 = numel(gp.nodes.functions.afid_areq0);

function gp = initcache(gp)
%INITCACHE Sets up fitness cache.
gp.fitness.cache = containers.Map('keytype','uint32','valuetype','any');
