function gp = gpcheck(gp)
%GPCHECK Perform pre-run error checks.
%
%   GP = GPCHECK(GP)
%
%   This function configures additional fields in the GP structure and
%   performs some error checking prior to a run.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPDEFAULTS, GPINIT

%set log zero warning state to off for duration of run and store current
%state so that it can be set back after run
gp.info.log0state = warning('query','MATLAB:log:logOfZero');
warning('off','MATLAB:log:logOfZero')

%automatically set number of inputs & number of train points if using
%multigene regression
if strncmpi(func2str(gp.fitness.fitfun),'regressmulti',12) || strncmpi(func2str(gp.fitness.fitfun),'binaryclass',11)
    [gp.userdata.numytrain, gp.nodes.inputs.num_inp] = size(gp.userdata.xtrain);
end

%some error checks
gp.runcontrol.pop_size = ceil(gp.runcontrol.pop_size);
if gp.runcontrol.pop_size < 2
    error('Population size must be at least 2. Set gp.runcontrol.pop_size >= 2')
end

gp.runcontrol.num_gen = ceil(gp.runcontrol.num_gen);
if gp.runcontrol.num_gen < 1
    error('GPTIPS must run for at least one generation. Set gp.runcontrol.num_gen >= 1');
end

gp.runcontrol.runs = ceil(gp.runcontrol.runs);
if gp.runcontrol.runs < 1
    error('GPTIPS must have at least 1 independent run. Set gp.runcontrol.runs >= 1');
end

gp.runcontrol.verbose = ceil(gp.runcontrol.verbose);
if gp.runcontrol.runs < 0
    error('The verbose parameter cannot be negative. Set gp.runcontrol.verbose >= 0');
end

gp.treedef.max_depth = ceil(gp.treedef.max_depth);
if gp.treedef.max_depth < 1
    error('gp.treedef.max_depth must be >= 1');
end

gp.treedef.max_mutate_depth = ceil(gp.treedef.max_mutate_depth);
if gp.treedef.max_mutate_depth < 1
    error('gp.treedef.max_mutate_depth must be >= 1');
end

if gp.treedef.max_mutate_depth > gp.treedef.max_depth
    gp.treedef.max_mutate_depth = gp.treedef.max_depth;
end

if ischar(gp.fitness.fitfun)
    gp.fitness.fitfun = str2func(gp.fitness.fitfun);
end

%for multigene regression check data
if strncmpi(func2str(gp.fitness.fitfun),'regressmulti',12) || strncmpi(func2str(gp.fitness.fitfun),'binaryclass',11)
    
    if ~isfield(gp.userdata,'xtrain') || ~isfield(gp.userdata,'ytrain') ...
            || isempty(gp.userdata.xtrain) || isempty(gp.userdata.ytrain)
        error('The fields gp.userdata.xtrain and gp.userdata.ytrain must be populated with data.');
    end
    
    if size(gp.userdata.xtrain,1) ~= size(gp.userdata.ytrain,1)
        error('There must be the same number of rows in gp.userdata.xtrain and gp.userdata.ytrain');
    end
    
    if any(~isfinite(gp.userdata.xtrain), 'all') || any(~isfinite(gp.userdata.ytrain), 'all')
        error('Non-finite values detected in gp.userdata.xtrain or gp.userdata.ytrain');
    end
    
    if any(~isreal(gp.userdata.xtrain), 'all') || any(~isreal(gp.userdata.ytrain), 'all')
        error('Complex values detected in gp.userdata.xtrain or gp.userdata.ytrain');
    end
    
end

gp.treedef.max_nodes = ceil(gp.treedef.max_nodes);
if gp.treedef.max_nodes < 1
    error('gp.treedef.max_nodes must be >= 1');
end

gp.treedef.build_method = ceil(gp.treedef.build_method);
if gp.treedef.build_method < 1 || gp.treedef.build_method > 3
    error('gp.treedef.build_method must be set as one of the following: 1 = full 2 = grow 3 = ramped half and half.');
end

if gp.selection.elite_fraction < 0 || gp.selection.elite_fraction > 1
    error('gp.selection.elite_fraction must be >= 0 and <= 1');
end

if sum([gp.operators.mutation.p_mutate gp.operators.crossover.p_cross gp.operators.directrepro.p_direct]) ~= 1
    error('gp crossover, mutation and direct reproduction probabilities must add up to 1');
end

if size(gp.operators.mutation.mutate_par,2) ~= 6 || ~isrow(gp.operators.mutation.mutate_par)
    error('gp.operators.mutation.mutate_par must be a row vector of length 6');
end

if sum(gp.operators.mutation.mutate_par) ~= 1
    error('The probabilities in gp.operators.mutation.mutate_par must add up to 1');
end

if gp.fitness.minimisation
    gp.fitness.minimisation = true;
end

if gp.fitness.complexityMeasure ~= 0 && gp.fitness.complexityMeasure ~= 1
    error('gp.fitness.complexityMeasure must be either 0 or 1');
end

if gp.nodes.const.p_ERC < 0 || gp.nodes.const.p_ERC > 1
    error('gp.nodes.const.p_ERC must be >= 0 and <= 1');
end

if gp.nodes.const.p_int < 0 || gp.nodes.const.p_int > 1
    error('gp.nodes.const.p_int must be >= 0 and <= 1');
end

gp.runcontrol.savefreq = ceil(gp.runcontrol.savefreq);
if gp.runcontrol.savefreq < 0
    error('gp.runcontrol.savefreq must be >= 0');
end

gp.nodes.inputs.num_inp = ceil(gp.nodes.inputs.num_inp);
if gp.nodes.inputs.num_inp < 0
    error('gp.nodes.inputs.num_inp must be >= 0');
end

gp.nodes.const.num_dec_places = ceil(gp.nodes.const.num_dec_places);

if gp.nodes.const.num_dec_places <= 0
    error('gp.nodes.const.num_dec_places must be > 0');
end

if any(size(gp.nodes.const.range) ~= [1 2])
    error('gp.nodes.const.range must be a 2 element row vector');
end

if gp.nodes.const.range(1) > gp.nodes.const.range(2)
    error('Invalid ERC range set in gp.nodes.const.range');
end

gp.nodes.const.rangesize = gp.nodes.const.range(2) - gp.nodes.const.range(1);

%if no inputs are defined then the probability of using an ERC must be set
%to 1 (if it is to be used at all).
if  gp.nodes.inputs.num_inp == 0 && gp.nodes.const.p_ERC > 0
    gp.nodes.const.p_ERC = 1;
end

%check that function node filenames have been supplied
if ~isfield(gp.nodes.functions,'name') || isempty(gp.nodes.functions.name) || ~iscell(gp.nodes.functions.name)
    error('Filenames for function nodes must be defined as a cell array in gp.nodes.functions.name ');
end


%loop through function nodes and check each is valid
for i = 1:length(gp.nodes.functions.name)
    
    status = exist(gp.nodes.functions.name{i});
    if (status ~= 2) && (status ~=3 ) && (status ~=5) && (status ~=6)
        error(['The function node "' gp.nodes.functions.name{i} ...
            '" cannot be found on the current path. Check your config file.']);
    end
    
end

%generate a function handle to the user defined fitness function
checkfile = exist(func2str(gp.fitness.fitfun),'file');
if ~(checkfile == 2 ) &&  ~(checkfile == 3) && ~(checkfile == 5) && ~(checkfile == 6)
    error('Cannot find a fitness function. In the config file this must be specified as a function handle in the field gp.fitness.fitfun');
end

%constant format control string
gp.nodes.const.format = [ '%0.' sprintf('%d',gp.nodes.const.num_dec_places) 'f' ];

if gp.genes.max_genes < 1
    error('gp.genes.max_genes must be > 0');
end

if gp.runcontrol.timeout <= 0
    error('The parameter gp.runcontrol.timeout must be >= 0');
end

if gp.selection.tournament.p_pareto < 0 || gp.selection.tournament.p_pareto > 1
    error('The parameter gp.selection.tournament.p_pareto must be >= 0 and <= 1');
end

if numel(gp.nodes.inputs.names) > gp.nodes.inputs.num_inp
    error('The supplied input variables ''name'' vector in gp.nodes.inputs.names contains more entries than input variables.');
end

%check and process user supplied variable aliases. For symbolic math
%toolbox lookup: remove whitespace, check entries are unique & remove
%markup
nameInds = find(~cellfun('isempty', gp.nodes.inputs.names));
numDefined = numel(nameInds);
numUnique = numel(unique(gp.nodes.inputs.names(nameInds)));

if numDefined ~= numUnique
    error('Entries of gp.nodes.inputs.names must be unique.');
end

%do some regex processing on names to get 'plain' versions
gp.nodes.inputs.namesPlain = gp.nodes.inputs.names;
gp.nodes.output.namePlain = gp.nodes.output.name;
gp.nodes.inputs.namesPlain(nameInds) = regexprep(gp.nodes.inputs.namesPlain(nameInds),'<.*?>','');
gp.nodes.inputs.namesPlain(nameInds) =  regexprep(gp.nodes.inputs.namesPlain(nameInds),'\s+','');
gp.nodes.output.namePlain =  regexprep(gp.nodes.output.namePlain,'<.*?>','');
gp.nodes.output.namePlain =  regexprep(gp.nodes.output.namePlain,'\s+','');
